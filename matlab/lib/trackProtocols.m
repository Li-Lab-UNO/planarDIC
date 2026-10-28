function trackProtocols(testRoot, settings)
% Track the deformations settings.track_deformations (Def01-Def19) against one common reference, the
% first frame of the reference deformation. The grid covers the region between the rest-state tine-tip
% lines; images are resampled to settings.target_px_per_mm, and the subset and grid step are set in mm.
% Each deformation is first matched from the reference frame to its own first frame (the bridge), then
% tracked through its frames; with settings.loading_only, only up to the end of the loading branch
% (util.loadingBranch). Writes DIC_results/DefNN_<name>.mat (positions in native image px).

cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
sc = util.readJSON(fullfile(testRoot, 'scale_calibration.json'));
f = min(1, settings.target_px_per_mm / sc.px_per_mm);           % resampling factor (downsampling only)
ppm = sc.px_per_mm * f;                                          % px/mm of the images tracked
radius = floor(round((settings.subset_mm * ppm - 1) / 2 + 0.5, 6));
step = floor(round(settings.step_mm * ppm + 0.5, 6));
toRes = @(v) (v + 0.5) * f - 0.5;                                % native px -> resampled px
toNat = @(v) (v + 0.5) / f - 0.5;
roi = latticeROI(toRes(cfg.tine_tips_px.z), toRes(cfg.tine_tips_px.theta), step);
[gx, gy] = meshgrid(roi(1):step:roi(3), roi(2):step:roi(4));
gx = gx'; gy = gy';                                              % points row by row (x fastest), as in Python
X0 = [gx(:), gy(:)]; shape = fliplr(size(gx));                   % shape = [ny nx]
maskSat = isfield(cfg, 'dic') && isfield(cfg.dic, 'mask_saturated') && cfg.dic.mask_saturated;
excl = zeros(0, 4);
if isfield(cfg, 'dic') && isfield(cfg.dic, 'exclude_px') && ~isempty(cfg.dic.exclude_px)
    excl = reshape(cfg.dic.exclude_px, [], 4);
end
prm = dicCore.trackParams(radius, maskSat, f);
prm.progress_s = settings.progress_s;
refFolder = util.defFolder(testRoot, settings.reference_deformation);
[tRef, pRef] = util.findFrames(fullfile(testRoot, 'Images', refFolder));
[img, use] = dicCore.loadFrame(pRef{1}, prm, toRes(excl));
ref = dicCore.buildReference(X0, radius, img, use);
util.log('  reference %s t = %.1f s; grid %d x %d points; subset %d px (%.2f mm), step %d px (%.2f mm); images x%.3f to %.1f px/mm%s', ...
         refFolder, tRef(1), shape(2), shape(1), 2 * radius + 1, (2 * radius + 1) / ppm, step, step / ppm, f, ppm, ...
         util.ternary(maskSat, '; saturated glints masked', ''));

P = util.phaseTiming(testRoot);
D = util.dataTable(testRoot);
defs = settings.track_deformations;
jobs = cell(numel(defs), 1);
for i = 1:numel(defs)
    fd = util.defFolder(testRoot, defs(i));
    [times, paths] = util.findFrames(fullfile(testRoot, 'Images', fd));
    branch = 'full';
    if settings.loading_only
        [~, t1] = util.loadingBranch(testRoot, defs(i), P, D);
        keep = times <= t1 + 0.05;
        times = times(keep); paths = paths(keep); branch = 'loading';
    end
    jobs{i} = struct('folder', fd, 'times', times, 'paths', {paths}, 'branch', branch);
end
outDir = fullfile(testRoot, 'DIC_results');
if ~isfolder(outDir), mkdir(outDir); end
meta = struct('X0', toNat(X0), 'shape', shape, 'step_px', step / f, 'radius', radius, 'resample', f, ...
    'px_per_mm', sc.px_per_mm, 'track_px_per_mm', ppm, 'subset_mm', (2 * radius + 1) / ppm, 'step_mm', step / ppm, ...
    'reference', sprintf('%s t=%.1fs', refFolder, tRef(1)), 'mask_saturated', maskSat, 'exclude_px', excl);
tAll = tic;
exclRes = toRes(excl);
[usePar, nw, nt] = parallelPlan(settings, numel(jobs));
if usePar
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= nw
        delete(pool); parpool(nw);
    end
    q = parallel.pool.DataQueue;                 % progress lines from the workers, printed as they come
    afterEach(q, @(m) util.log('%s', m));
    report = @(m, varargin) send(q, m);         % one line per update: no overwriting with several workers
    prm.progress_s = max(prm.progress_s, 30);
    util.log('  tracking %d deformations, %d at a time, %d thread(s) each', numel(jobs), nw, nt);
    parfor i = 1:numel(jobs)
        maxNumCompThreads(nt);                   % a worker uses one core unless told otherwise
        trackOne(jobs{i}, ref, X0, prm, exclRes, meta, outDir, toNat, report);
    end
else
    for i = 1:numel(jobs)
        util.log('  [%d/%d] %s: %d frames (%s)', i, numel(jobs), jobs{i}.folder, numel(jobs{i}.paths), jobs{i}.branch);
        trackOne(jobs{i}, ref, X0, prm, exclRes, meta, outDir, toNat, @util.progress);
    end
end
util.log('  all %d deformations tracked in %.1f min', numel(jobs), toc(tAll) / 60);
end


function trackOne(job, ref, X0, prm, exclRes, meta, outDir, toNat, report)
% Bridge from the reference frame to this deformation's first frame, then track its frames.
img0 = dicCore.loadFrame(job.paths{1}, prm);
[p, c, cv] = dicCore.icgn(ref, img0, zeros(size(X0, 1), 6), prm.max_iter, prm.tol);
good = cv & c >= prm.zncc_min;
p(~isfinite(p)) = 0;
p = dicCore.fillFromNeighbors(X0, p, good);
Xstart = X0 + p(:, [1 4]);
prm.exclude_px = dicCore.shiftRects(exclRes, X0, Xstart, good);   % excluded regions into this first frame
R = dicCore.trackSequence(job.paths, Xstart, prm, job.times, strtok(job.folder, '_'), report);
R.valid = R.valid & good';
track = meta;
track.times = job.times(:);
[~, names] = cellfun(@fileparts, job.paths, 'UniformOutput', false);
track.frames = names(:);
track.posx = toNat(R.posx); track.posy = toNat(R.posy);
track.valid = R.valid; track.zncc = R.zncc; track.keyframes = R.keyframes;
track.bridge_zncc = c; track.branch = job.branch;
save(fullfile(outDir, [job.folder '.mat']), '-struct', 'track');
report(sprintf('    %s done: valid per frame median %.1f%%, min %.1f%%; bridge valid %.1f%%', job.folder, ...
               100 * median(mean(R.valid, 2)), 100 * min(mean(R.valid, 2)), 100 * mean(good)));
end


function roi = latticeROI(zTips, thTips, step)
% Grid extents [x0 y0 x1 y1]: a lattice of `step` centered between the tine-tip lines, >= 5 px inside them.
roi = zeros(1, 4);
[roi(1), roi(3)] = axisRange(zTips(1), zTips(2));
[roi(2), roi(4)] = axisRange(thTips(1), thTips(2));
    function [x0, x1] = axisRange(lo, hi)
        a = lo + 5; b = hi - 5;
        n = floor((b - a) / step);
        x0 = round(a + ((b - a) - n * step) / 2);
        x1 = x0 + n * step;
    end
end


function [usePar, nw, nt] = parallelPlan(settings, nJobs)
% Whether to track several deformations at once (settings.use_parallel: 'auto' does with the Parallel
% Computing Toolbox and at least 16 cores), how many (settings.parallel_workers: 'auto' is half the cores,
% at most one per deformation) and with how many threads each (settings.parallel_threads: 'auto' shares
% all the processor's threads, hyperthreads included, among the workers).
cores = feature('numcores');
mode = settings.use_parallel;
if ischar(mode) || isstring(mode), usePar = cores >= 16; else, usePar = logical(mode); end
usePar = usePar && nJobs > 1 && canUseParallel();
nw = settings.parallel_workers;
if ischar(nw) || isstring(nw), nw = max(2, floor(cores / 2)); end
nw = max(1, min(nw, nJobs));
nt = settings.parallel_threads;
if ischar(nt) || isstring(nt), nt = max(1, floor(logicalProcessors(cores) / nw)); end
end


function n = logicalProcessors(cores)
% The processor's threads (2 per core with hyperthreading); the cores if that can't be read.
n = cores;
try
    n = max(cores, java.lang.Runtime.getRuntime().availableProcessors());
catch
end
end


function ok = canUseParallel()
% True when the Parallel Computing Toolbox is installed and licensed.
ok = false;
try
    ok = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
catch
end
end
