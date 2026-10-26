function changed = artifactCheck(testRoot, settings, reportOnly)
% Find local tracking errors and exclude their cause from the matching. Some image features move
% differently from the tissue: a speckle particle that slides over it, or a saturated glare patch whose
% edges shift with load. Subsets overlapping one follow it. They are found and confirmed automatically:
%   1. Candidates: at the peak of the reference deformation (Def11), each grid point's displacement is
%      compared with the median of its 9 x 9 neighborhood. Clusters of at least 6 connected points above
%      max(60 um, 8 x the median residual), peaking at >= 100 um and >= 1 mm from every tine-tip line
%      (tine contact zones have steep real gradients), are candidates.
%   2. Confirmation: small textured patches of the rest frame are matched into the peak frame by plain
%      normalized cross-correlation, independently of the tracker. A candidate is confirmed when the
%      tracker disagrees with these matches by >= 3 px (median) inside it and >= 3 times more than in a
%      ring 1-2 mm around it.
%   3. Fix: each confirmed cluster's bounding box (reference-frame px, grown by half a grid step) is added
%      to dic.exclude_px in test_config.json. Those pixels are left out of the matching; the test then has
%      to be tracked again (runDIC does this).
% Returns true when test_config.json was changed. Writes DIC_results/artifact_check.png.

if nargin < 3, reportOnly = false; end
C.RES_MIN_UM = 60; C.RES_FACTOR = 8; C.MIN_POINTS = 6; C.PEAK_UM = 100; C.EDGE_MM = 1;
C.CONFIRM_PX = 3; C.CONFIRM_RATIO = 3;

cfgFile = fullfile(testRoot, 'test_config.json');
cfg = util.readJSON(cfgFile);
sc = util.readJSON(fullfile(testRoot, 'scale_calibration.json'));
mpp = sc.mm_per_px;
folder = util.defFolder(testRoot, settings.reference_deformation);
R = load(fullfile(testRoot, 'DIC_results', [folder '.mat']));
[~, paths] = util.findFrames(fullfile(testRoot, 'Images', folder));
paths = paths(ismember(cellfun(@(p) fileName(p), paths, 'UniformOutput', false), R.frames));
ny = R.shape(1); nx = R.shape(2);
grid2 = @(v) reshape(v, nx, ny)';               % points are stored row by row (x fastest)

% peak frame: largest mean displacement of the valid points
ux = R.posx - R.X0(:, 1)'; uy = R.posy - R.X0(:, 2)';
mag = hypot(ux, uy); mag(~R.valid) = NaN;
[~, k] = max(mean(mag, 2, 'omitnan'));
% residual from the 9 x 9 neighborhood median (um)
res = zeros(ny, nx);
for comp = {ux, uy}
    U = comp{1}(k, :); U(~R.valid(k, :)) = NaN;
    U = grid2(U);
    res = res + (U - nanMedianFilter(U, 4)) .^ 2;
end
res = sqrt(res) * mpp * 1000;
tips = cfg.tine_tips_px;
X0x = grid2(R.X0(:, 1)); X0y = grid2(R.X0(:, 2));
z = (X0x - tips.z(1)) * mpp; th = (tips.theta(2) - X0y) * mpp;
Lz = (tips.z(2) - tips.z(1)) * mpp; Lt = (tips.theta(2) - tips.theta(1)) * mpp;
edge = min(min(z, Lz - z), min(th, Lt - th));
bg = median(res(:), 'omitnan');
thr = max(C.RES_MIN_UM, C.RES_FACTOR * bg);
r0 = res; r0(isnan(r0)) = 0;
[lab, n] = label4(r0 > thr);
g0 = imageOps.gaussBlur(imageOps.readGray(paths{1}), 1);
gk = imageOps.gaussBlur(imageOps.readGray(paths{k}), 1);
v = R.valid(k, :)';
u = [ux(k, :)', uy(k, :)'];
interpU = scatteredInterpolant(R.X0(v, 1), R.X0(v, 2), u(v, 1), 'linear', 'none');
interpV = scatteredInterpolant(R.X0(v, 1), R.X0(v, 2), u(v, 2), 'linear', 'none');
found = struct('z', {}, 'theta', {}, 'points', {}, 'peak_um', {}, 'confirmed', {}, 'why', {}, 'rect', {}, ...
               'match_in', {}, 'match_ring', {});
for i = 1:n
    m = lab == i;
    c = struct('z', mean(z(m)), 'theta', mean(th(m)), 'points', sum(m(:)), 'peak_um', max(res(m)), ...
               'confirmed', false, 'why', '', 'rect', [], 'match_in', NaN, 'match_ring', NaN);
    if c.points < C.MIN_POINTS || c.peak_um < C.PEAK_UM, continue; end      % small or weak: noise
    if min(edge(m)) < C.EDGE_MM
        c.why = sprintf('within %g mm of a tine line', C.EDGE_MM);
        found(end + 1) = c; %#ok<AGROW>
        continue
    end
    xs = X0x(m); ys = X0y(m);
    x0 = min(xs) - R.step_px / 2; y0 = min(ys) - R.step_px / 2;
    x1 = max(xs) + R.step_px / 2; y1 = max(ys) + R.step_px / 2;
    [PX, PY] = meshgrid(x0:15:x1 - 1e-9, y0:15:y1 - 1e-9);
    [rr, aa] = meshgrid((1:0.25:2) / mpp, (0:23) * 2 * pi / 24);
    ring = [mean(xs) + rr(:) .* cos(aa(:)), mean(ys) + rr(:) .* sin(aa(:))];
    dIn = patchDisagreement(g0, gk, [PX(:), PY(:)], interpU, interpV);
    dRing = patchDisagreement(g0, gk, ring, interpU, interpV);
    c.match_in = median(dIn); c.match_ring = median(dRing);
    c.rect = [floor(x0), floor(y0), ceil(x1), ceil(y1)];
    if numel(dIn) < 3 || numel(dRing) < 5
        c.why = 'too few textured patches to confirm';
    elseif c.match_in >= C.CONFIRM_PX && c.match_in >= C.CONFIRM_RATIO * max(c.match_ring, 0.5)
        c.confirmed = true; c.why = 'tracker disagrees with independent patch matching';
    else
        c.why = 'patch matching agrees with the tracker';
    end
    found(end + 1) = c; %#ok<AGROW>
end

util.log('  %s peak frame %d: median residual %.1f um, threshold %.0f um', folder, k, bg, thr);
for c = found
    extra = '';
    if isfinite(c.match_in), extra = sprintf('; patch matching vs tracker %.1f px inside, %.1f px around', c.match_in, c.match_ring); end
    util.log('  (%.1f, %.1f) mm: %d points, peak %.0f um: %s, %s%s', c.z, c.theta, c.points, c.peak_um, ...
             util.ternary(c.confirmed, 'CONFIRMED', 'not confirmed'), c.why, extra);
end
if isempty(found), util.log('  no tracking artifacts found'); end
excl = zeros(0, 4);
if isfield(cfg, 'dic') && isfield(cfg.dic, 'exclude_px') && ~isempty(cfg.dic.exclude_px), excl = reshape(cfg.dic.exclude_px, [], 4); end
plotCheck(fullfile(testRoot, 'DIC_results', 'artifact_check.png'), z, th, res, thr, found, cfg.sample, excl, tips, mpp);

changed = false;
conf = found([found.confirmed]);
if isempty(conf) || reportOnly, return; end
rects = zeros(0, 4);
if isfield(cfg.dic, 'exclude_px') && ~isempty(cfg.dic.exclude_px), rects = reshape(cfg.dic.exclude_px, [], 4); end
added = {};
for c = conf
    r = c.rect;
    overlap = any(min(r(3), rects(:, 3)) > max(r(1), rects(:, 1)) & min(r(4), rects(:, 4)) > max(r(2), rects(:, 2)));
    if ~overlap
        rects(end + 1, :) = r; %#ok<AGROW>
        added{end + 1} = sprintf('[%d, %d, %d, %d] at (%.1f, %.1f) mm', r, c.z, c.theta); %#ok<AGROW>
    end
end
if isempty(added)
    util.log('  the confirmed artifacts already lie in excluded regions: nothing to change');
    return
end
cfg.dic.exclude_px = num2cell(rects, 2);        % a list of [x0 y0 x1 y1] rows, as Python writes it
note = sprintf('artifact check excluded %s (tracker disagreed with independent patch matching)', strjoin(added, '; '));
if isfield(cfg.dic, 'note') && ~isempty(cfg.dic.note), cfg.dic.note = [cfg.dic.note '; ' note]; else, cfg.dic.note = note; end
util.writeJSON(cfgFile, cfg);
util.log('  excluded %s: the test is tracked again', strjoin(added, '; '));
changed = true;
end


function d = patchDisagreement(g0, gk, pts, interpU, interpV)
% For textured 25 px patches of the rest frame g0 centered at pts (px): distance (px) between where plain
% normalized cross-correlation finds them in gk and where the tracked displacement puts them.
S = 45; [h, w] = size(g0);
d = [];
for q = 1:size(pts, 1)
    x = round(pts(q, 1)); y = round(pts(q, 2));
    if x - 12 < 0 || y - 12 < 0 || x + 12 > w - 1 || y + 12 > h - 1, continue; end
    tpl = g0(y - 11:y + 13, x - 11:x + 13);
    if std(tpl(:), 1) < 10, continue; end
    ug = [interpU(x, y), interpV(x, y)];
    if ~all(isfinite(ug)), continue; end
    gx = fix(x + ug(1)); gy = fix(y + ug(2));
    if gx - 12 - S < 0 || gy - 12 - S < 0 || gx + 12 + S > w - 1 || gy + 12 + S > h - 1, continue; end
    win = gk(gy - 11 - S:gy + 13 + S, gx - 11 - S:gx + 13 + S);
    [mx, lx, ly] = matchTemplate(win, tpl);
    if mx < 0.8, continue; end
    d(end + 1) = hypot(gx - S + lx - (x + ug(1)), gy - S + ly - (y + ug(2))); %#ok<AGROW>
end
end


function [best, lx, ly] = matchTemplate(win, tpl)
% Normalized cross-correlation (correlation coefficient) of template tpl at every position inside win;
% returns the best value and its 0-based offset (x, y) in win.
t = tpl - mean(tpl(:));
nt = numel(t);
num = conv2(win, rot90(t, 2), 'valid');
s1 = conv2(win, ones(size(t)), 'valid');
s2 = conv2(win .^ 2, ones(size(t)), 'valid');
den = sqrt(sum(t(:) .^ 2) * max(s2 - s1 .^ 2 / nt, 0));
r = num ./ den; r(den == 0) = -1;
[best, i] = max(r(:));
[iy, ix] = ind2sub(size(r), i);
lx = ix - 1; ly = iy - 1;
end


function M = nanMedianFilter(U, r)
% Median of the finite values in the (2r+1) x (2r+1) neighborhood of every element (edges replicated).
[ny, nx] = size(U);
iy = min(max((1 - r:ny + r), 1), ny); ix = min(max((1 - r:nx + r), 1), nx);
P = U(iy, ix);
stack = zeros(ny, nx, (2 * r + 1) ^ 2);
q = 0;
for dy = 0:2 * r
    for dx = 0:2 * r
        q = q + 1;
        stack(:, :, q) = P(1 + dy:ny + dy, 1 + dx:nx + dx);
    end
end
M = median(stack, 3, 'omitnan');
end


function [lab, n] = label4(mask)
% Connected components of a logical matrix (4-connected), numbered 1..n.
lab = zeros(size(mask)); n = 0;
[ny, nx] = size(mask);
for start = find(mask)'
    if lab(start), continue; end
    n = n + 1; lab(start) = n; stack = start;
    while ~isempty(stack)
        p = stack(end); stack(end) = [];
        [iy, ix] = ind2sub([ny nx], p);
        nb = [iy - 1 ix; iy + 1 ix; iy ix - 1; iy ix + 1];
        nb = nb(nb(:, 1) >= 1 & nb(:, 1) <= ny & nb(:, 2) >= 1 & nb(:, 2) <= nx, :);
        for j = 1:size(nb, 1)
            pj = sub2ind([ny nx], nb(j, 1), nb(j, 2));
            if mask(pj) && ~lab(pj), lab(pj) = n; stack(end + 1) = pj; end %#ok<AGROW>
        end
    end
end
end


function plotCheck(path, z, th, res, thr, found, sample, excl, tips, mpp)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 760]);
pcolor(z, th, res); shading flat; axis equal tight; hold on
colormap(hot); clim([0 max(2 * thr, 100)]);
cb = colorbar; cb.Label.String = '|displacement - neighborhood median| at the peak (\mum)';
green = [0.18 0.8 0.44];
for q = 1:size(excl, 1)                         % excluded regions (not tracked: blank in the map)
    r = excl(q, :);
    zr = ([r(1) r(3)] - tips.z(1)) * mpp; tr = (tips.theta(2) - [r(2) r(4)]) * mpp;
    plot(zr([1 2 2 1 1]), tr([1 1 2 2 1]), '--', 'Color', green, 'LineWidth', 1.5);
    text(mean(zr), min(tr) - 0.15, 'excluded', 'Color', green, 'FontSize', 9, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', 'FontWeight', 'bold');
end
zMid = mean(z(:), 'omitnan');
for c = found
    col = green; if ~c.confirmed, col = [0.62 0.77 0.96]; end
    plot(c.z, c.theta, 'o', 'MarkerSize', 20, 'LineWidth', 2, 'Color', col);
    txt = sprintf('%d pts, %.0f \\mum\n%s', c.points, c.peak_um, c.why);
    if c.z < zMid
        text(c.z + 0.35, c.theta, txt, 'Color', col, 'FontSize', 8);
    else
        text(c.z - 0.35, c.theta, txt, 'Color', col, 'FontSize', 8, 'HorizontalAlignment', 'right');
    end
end
xlabel('z, axial (mm)'); ylabel('\theta, circumferential (mm)');
title({strrep(sample, '_', ' '), 'tracking artifact check: green = excluded artifact, blue = checked, not an artifact'}, 'FontWeight', 'normal');
print(fig, path, '-dpng', '-r100');
close(fig);
end


function n = fileName(p)
[~, n] = fileparts(p);
end
