function cfg = findTineTips(testRoot, label, settings)
% Find the 20 tine tips in the rest frame (first image) of the reference deformation (Def11) and write the
% tine settings into test_config.json: tine_tips_px (the rest-state tip line of each rake), the tips of
% every tine (calibration_tines_px), the tracking regions, the sample label and the glint check.
%
% Rods: each rake's five rods are long straight dark lines in the background band between the image edge
% and the tissue. In every column (left/right rakes) or row (top/bottom rakes) of that band the profile
% across the rods is searched for narrow dark dips; dips that line up across the band are fitted with a
% straight line per rod. Tips: each rod is followed inward along its line; the tip is the end of the
% rod-dark stretch at the largest sustained step up in brightness (rod to tissue). Tips far from their
% rake's median are searched again. Every tip is then refined with the scale calibration's tip finder.
% A self-check (five rods per rake, tips of a rake within 35 px of each other, a clear step at every tip)
% stops the run if it fails: check DIC_results/config_tine_tips.png.

C.N_RODS = 5;
C.BANDS = [0.10 0.27; 0.04 0.18; 0.08 0.20; 0.16 0.30];
C.SIGMAS = [4 8];
C.MAX_SHIFT = 6; C.RAKE_TOL = 20; C.SETTLE_PX = 3; C.MAX_GAP = 12; C.SPREAD_TOL = 35; C.STEP_PX = 30;
C.RETRY = [0.3 20; 0.5 12; 0.5 20];
C.SAT_LEVEL = 240; C.SAT_MAX = 10;
rakeNames = {'left', 'right', 'top', 'bottom'};
rakeHor = [true true false false];
rakeOut = [-1 1 -1 1];

refDef = settings.reference_deformation;
folder = util.defFolder(testRoot, refDef);
[~, paths] = util.findFrames(fullfile(testRoot, 'Images', folder));
gray8 = imageOps.readGray(paths{1});
g = imageOps.gaussBlur(gray8, 1.0);
[H, W] = size(g);

% ---- detect rods and tips (or, with settings.tine_tips = 'click', ask for all of them)
rakes = struct(); problems = {};
detectRakes = 1:4;
clicked = false(1, 4);                           % rakes whose tips were clicked: kept as clicked
if strcmpi(settings.tine_tips, 'click')
    for r = 1:4
        rakes.(rakeNames{r}) = clickRake(gray8, rakeNames{r}, rakeHor(r), rakeOut(r));
    end
    detectRakes = []; clicked(:) = true;
end
for r = detectRakes
    hor = rakeHor(r); out = rakeOut(r);
    n = W * hor + H * ~hor;
    best = []; bestLines = [];
    for sig = C.SIGMAS
        for b = 1:size(C.BANDS, 1)
            f0 = C.BANDS(b, 1); f1 = C.BANDS(b, 2);
            if out < 0, bd = [fix(f0 * n), fix(f1 * n)]; else, bd = [fix((1 - f1) * n), fix((1 - f0) * n)]; end
            ls = findRods(g, hor, bd, sig, C);
            if isempty(best) || size(ls, 1) > size(bestLines, 1), best = bd; bestLines = ls; end
            if size(ls, 1) == C.N_RODS, break; end
        end
        if size(bestLines, 1) == C.N_RODS, break; end
    end
    if size(bestLines, 1) ~= C.N_RODS
        problems{end + 1} = sprintf('%s rake: found %d rods, expected %d', rakeNames{r}, size(bestLines, 1), C.N_RODS); %#ok<AGROW>
        continue
    end
    if out < 0, start = best(1); else, start = best(2); end
    found = zeros(C.N_RODS, 5);
    for k = 1:C.N_RODS
        [found(k, 1), found(k, 2), found(k, 3), found(k, 4), found(k, 5)] = ...
            followRod(g, bestLines(k, :), hor, out, start, fix(n / 2), 0.3, C.MAX_GAP, false, [], C);
    end
    med = median(found(:, 1));
    tips = zeros(C.N_RODS, 2);
    for k = 1:C.N_RODS
        f = found(k, :);
        if abs(f(1) - med) > C.RAKE_TOL        % tines are the same length: search again near the median
            alts = [];
            for q = 1:size(C.RETRY, 1)
                for lastStep = [false true]
                    [a1, a2, a3, a4, a5] = followRod(g, bestLines(k, :), hor, out, med + out * 80, med - out * 80, ...
                        C.RETRY(q, 1), C.RETRY(q, 2), lastStep, f(4:5), C);
                    alts(end + 1, :) = [a1 a2 a3 a4 a5]; %#ok<AGROW>
                end
            end
            [~, j] = min(abs(alts(:, 1) - med));
            f = alts(j, :);
        end
        t = f(1); c = f(2); step = f(3); rodGray = f(4); bg = f(5);
        weak = step < max(25, 0.4 * (bg - rodGray));
        if abs(t - med) > C.RAKE_TOL || weak     % last resort: the calibration's tip finder from the median
            cands = [];
            for dk = [-20 -10 0 10 20]
                [tn, con] = tipFrom(g, c, med + dk, hor, out);
                if isfinite(tn) && con > 30 && abs(tn - med) <= C.SPREAD_TOL, cands(end + 1, :) = [tn con]; end %#ok<AGROW>
            end
            if ~isempty(cands)
                [~, j] = min(abs(cands(:, 1) - med));
                t = cands(j, 1); weak = false;
            end
        end
        if weak
            problems{end + 1} = sprintf('%s rake, rod at %.0f px: weak tip step (%.0f gray levels)', rakeNames{r}, c, step); %#ok<AGROW>
        end
        tips(k, :) = [round(c), round(t)];
    end
    spread = max(abs(tips(:, 2) - median(tips(:, 2))));
    if spread > C.SPREAD_TOL
        problems{end + 1} = sprintf('%s rake: tips spread %.0f px from their median (tines are the same length)', rakeNames{r}, spread); %#ok<AGROW>
    end
    rakes.(rakeNames{r}) = struct('horizontal', hor, 'out', out, 'tips', tips);
end
if ~isempty(problems)
    util.log('  Tine-tip self-check failed:\n    %s', strjoin(problems, '\n    '));
    failed = unique(cellfun(@(q) strtok(q), problems, 'UniformOutput', false));
    cfgFile = fullfile(testRoot, 'test_config.json');
    old = struct();
    if isfile(cfgFile)
        c0 = util.readJSON(cfgFile);
        if isfield(c0, 'calibration_tines_px'), old = c0.calibration_tines_px; end
    end
    for r = find(ismember(rakeNames, failed))
        if isfield(old, rakeNames{r})                % tines already in the config: start from them
            rakes.(rakeNames{r}) = old.(rakeNames{r});
            util.log('  %s rake: using the tines already in test_config.json', rakeNames{r});
        elseif usejava('desktop')                   % ask for the tips
            near = [];
            if isfield(rakes, rakeNames{r}), near = rakes.(rakeNames{r}).tips; end
            rakes.(rakeNames{r}) = clickRake(gray8, rakeNames{r}, rakeHor(r), rakeOut(r), near);
            clicked(r) = true;
            util.log('  %s rake: using the tips clicked', rakeNames{r});
        else
            overlay(gray8, rakes, [], fullfile(testRoot, 'DIC_results', 'config_tine_tips.png'));
            error('runDIC:tines', ['Tine-tip self-check failed for the %s rake. Run runDIC in the MATLAB ' ...
                  'desktop to click its tips, or see DIC_results/config_tine_tips.png.'], rakeNames{r});
        end
    end
end

% ---- refine every tip with the calibration's tip finder (tip coordinate only)
g08 = imageOps.gaussBlur(gray8, 0.8);
[rakes, keptBy] = refineRakes(g08, rakes, rakeNames(~clicked), C, struct());
for r = find(clicked), keptBy.(rakeNames{r}) = C.N_RODS; end   % clicked tips are kept as clicked
tipsPx = tipLines(rakes);

% ---- ask whether the tips are right, and let the wrong rakes be clicked (desktop only)
if settings.review_tines && usejava('desktop')
    saved = hasTips(testRoot);                         % tips saved before: Skip keeps them
    while true
        answer = askTips(gray8, rakes, tipsPx, rakeNames, saved);
        if strcmp(answer, 'Cancel'), util.cancel('the tine tips'); end
        if strcmp(answer, 'Skip')
            util.log('Tine tips: skipped, the tips saved before are kept');
            cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
            return
        end
        r = find(strcmp(answer, clickLabels(rakeNames)));
        if isempty(r), break; end                  % Continue
        rakes.(rakeNames{r}) = clickRake(gray8, rakeNames{r}, rakeHor(r), rakeOut(r), rakes.(rakeNames{r}).tips);
        keptBy.(rakeNames{r}) = C.N_RODS;              % kept as clicked
        tipsPx = tipLines(rakes);
    end
end
z0 = tipsPx.z(1); z1 = tipsPx.z(2); y0 = tipsPx.theta(1); y1 = tipsPx.theta(2);
region = gray8(y0 + 1:y1, z0 + 1:z1);
sat = 100 * mean(region(:) >= C.SAT_LEVEL);

% ---- write test_config.json (other settings already in it are kept)
cfgFile = fullfile(testRoot, 'test_config.json');
if isfile(cfgFile), cfg = util.readJSON(cfgFile); else, cfg = struct(); end
if ~isempty(label) || ~isfield(cfg, 'sample'), cfg.sample = util.sampleName(testRoot, label); end
cfg.roi_px = struct('standard', [z0 + 40, y0 + 40, z1 - 40, y1 - 40], 'tine_region', [z0 + 10, y0 + 10, z1 - 10, y1 - 10]);
cfg.tine_tips_px = tipsPx;
cfg.calibration_tines_px = rakes;
if ~isfield(cfg, 'reference_deformation'), cfg.reference_deformation = sprintf('Def%02d', refDef); end
if ~isfield(cfg, 'dic'), cfg.dic = struct(); end
if ~isfield(cfg.dic, 'mask_saturated')
    cfg.dic.mask_saturated = sat > 0 && sat < C.SAT_MAX;
    if sat == 0, why = 'no saturated pixels'; elseif sat < C.SAT_MAX, why = 'glints masked';
    else, why = 'overexposed: masking would remove too much of the speckle'; end
    cfg.dic.note = sprintf('rest frame: %.3f%% of the tine region at gray >= %d (%s)', sat, C.SAT_LEVEL, why);
end
util.writeJSON(cfgFile, cfg);
overlay(gray8, rakes, tipsPx, fullfile(testRoot, 'DIC_results', 'config_tine_tips.png'));
util.log('Tine tips: z %d-%d px, theta %d-%d px (%d of 20 refined); saturated %.2f%% of the tine region', ...
         z0, z1, y0, y1, 20 - sum(cellfun(@(f) keptBy.(f), fieldnames(keptBy))), sat);
end


% ======================================================================================================
function lines = findRods(g, hor, band, sigma, C)
% Straight-line fits [a b score] (c = a + b s, s along the rods) of the rake's rods inside the band.
[H, W] = size(g);
if hor, nAcross = H; else, nAcross = W; end
lo = fix(0.2 * nAcross); hi = fix(0.8 * nAcross);
hits = [];
kb = fix(15 * sigma); if mod(kb, 2) == 0, kb = kb + 1; end
for s = band(1):6:band(2) - 1
    if hor, prof = mean(g(lo + 1:hi, s - 1:s + 3), 2)'; else, prof = mean(g(s - 1:s + 3, lo + 1:hi), 1); end
    prof = blur1(prof, sigma);
    base = box1(prof, kb);
    [pk, pr] = findPeaks(base - prof, 20, [4 40], 25);
    hits = [hits; repmat(s, numel(pk), 1), lo + pk(:) - 1, pr(:)]; %#ok<AGROW>
end
lines = zeros(0, 3);
if isempty(hits), return; end
used = false(size(hits, 1), 1);
nMin = (band(2) - band(1)) / 6;
for c0 = sort(hits(:, 2))'
    m = ~used & abs(hits(:, 2) - c0) < 10;
    if sum(m) < 0.5 * nMin, continue; end
    for it = 1:3
        pf = polyfit(hits(m, 1), hits(m, 2), 1);
        m = ~used & abs(hits(:, 2) - (pf(2) + pf(1) * hits(:, 1))) < 4;
    end
    if sum(m) >= 0.6 * nMin
        used = used | m;
        pf = polyfit(hits(m, 1), hits(m, 2), 1);
        lines(end + 1, :) = [pf(2), pf(1), sum(hits(m, 3))]; %#ok<AGROW>
    end
end
[~, o] = sort(-lines(:, 3)); lines = lines(o(1:min(C.N_RODS, end)), :);
[~, o] = sort(lines(:, 1) + lines(:, 2) * mean(band)); lines = lines(o, :);
end


function [tip, cAt, step, rodGray, bg] = followRod(g, line, hor, out, sStart, sStop, dark, maxGap, lastStep, levels, C)
% Follow the rod inward from sStart toward sStop (see findTineTips). Returns the tip (coordinate along the
% rod), the centerline at the tip, the brightness step, and the rod and background gray levels.
a = line(1); b = line(2);
off = -15:15; core = abs(off) <= 6; side = abs(off) >= 10;
sAll = sStart - out * (0:ceil(abs(sStop - sStart)));     % from sStart toward sStop, sStop excluded
if out < 0, sAll = sAll(sAll < sStop); else, sAll = sAll(sAll > sStop); end
shift = 0; cs = zeros(size(sAll)); cen = cs; sid = cs;
for q = 1:numel(sAll)
    s = sAll(q);
    c = a + b * s + shift;
    P = sampleLine(g, repmat(s, 1, numel(off)), c + off, hor);
    P = conv(P, ones(1, 3) / 3, 'same');
    Pc = P; Pc(abs(off) > 3) = inf;
    [~, i] = min(Pc);
    shift = min(max(shift + 0.3 * off(i), -C.MAX_SHIFT), C.MAX_SHIFT);
    cs(q) = c; cen(q) = min(P(core)); sid(q) = mean(P(side));
end
if isempty(levels)
    rodGray = median(cen(1:min(60, end))); bg = median(sid(1:min(60, end)));
else
    rodGray = levels(1); bg = levels(2);
end
contrast = max([bg - rodGray, 0.5 * (util.percentile(cen, 90) - rodGray), 30]);
S = C.STEP_PX;
k = (S + 5):(numel(cen) - S - 6);                 % 0-based q
score = -inf(size(k));
for j = 1:numel(k)
    q = k(j) + 1;                                 % 1-based
    before = median(cen(q - S - 3:q - 4));
    after = median(cen(q + 3:q + S + 2));
    if before < rodGray + 0.35 * contrast, score(j) = after - before; end
end
step = max(score);
if lastStep, j = find(score >= 0.6 * step, 1, 'last'); else, [~, j] = max(score); end
q = k(j);                                         % 0-based
on = cen < rodGray + dark * max(contrast, step);
last = max(q - 40, 0); gap = 0;
for i = max(q - 40, 0):min(q + 40, numel(cen)) - 1
    if on(i + 1), last = i; gap = 0;
    else
        gap = gap + 1;
        if gap > maxGap, break; end
    end
end
tip = sAll(last + 1); cAt = cs(last + 1);
end


function [tt, con] = tipFrom(g, c, t, hor, out)
% Tip found by the calibration's tip finder started at tip coordinate t on the rod at c (6 passes).
cc = c; tt = t; con = 0;
for it = 1:6
    [cc, tn, con] = calFindTip(g, cc, tt, hor, out);
    if ~isfinite(tn), tt = NaN; return; end
    tt = tn;
end
end


function v = sampleLine(g, along, across, hor)
if hor, v = imageOps.sampleLinear(g, along, across, 'replicate');
else, v = imageOps.sampleLinear(g, across, along, 'replicate'); end
end


function y = blur1(x, sigma)
% 1D Gaussian blur of a row vector, mirrored borders.
k = round(sigma * 8 + 1); if mod(k, 2) == 0, k = k + 1; end
r = (k - 1) / 2; w = exp(-(-r:r) .^ 2 / (2 * sigma ^ 2)); w = w / sum(w);
y = conv(x(imageOps.reflectIdx(numel(x), r)), w, 'valid');
end


function y = box1(x, k)
% 1D moving mean of width k (odd), mirrored borders.
r = (k - 1) / 2;
y = conv(x(imageOps.reflectIdx(numel(x), r)), ones(1, k) / k, 'valid');
end


function [rakes, keptBy] = refineRakes(g08, rakes, names, C, keptBy)
% Refine the tips of rakes `names` with the calibration's tip finder (4 passes; tip coordinate only, the
% centerline it finds sits on a rod edge). A refinement is kept when its last two passes agree within
% SETTLE_PX, the rod-to-tissue contrast exceeds 30 gray levels and it moves the tip by at most 25 px
% (centerline at most 10 px); otherwise the tip stays as found. Clicked tips are not refined.
% keptBy.(rake) = tips kept unrefined.
for r = 1:numel(names)
    R = rakes.(names{r});
    kept = 0;
    for k = 1:size(R.tips, 1)
        c0 = R.tips(k, 1); t0 = R.tips(k, 2);
        c = c0; t = t0; ts = []; con = 0;
        for it = 1:4
            [c, tn, con] = calFindTip(g08, c, t, R.horizontal, R.out);
            if isfinite(tn), t = tn; ts(end + 1) = tn; end %#ok<AGROW>
        end
        settled = numel(ts) >= 2 && abs(ts(end) - ts(end - 1)) <= C.SETTLE_PX;
        if settled && con > 30 && abs(t - t0) <= 25 && abs(c - c0) <= 10
            R.tips(k, :) = [c0, round(t)];
        else
            kept = kept + 1;
        end
    end
    rakes.(names{r}) = R;
    keptBy.(names{r}) = kept;
end
end


function tipsPx = tipLines(rakes)
% Rest-state tine-tip lines: the median tip of each rake.
medTip = @(s) round(median(rakes.(s).tips(:, 2)));
tipsPx = struct('z', [medTip('left'), medTip('right')], 'theta', [medTip('top'), medTip('bottom')]);
end


function R = clickRake(gray8, name, hor, out, near)
% Ask for the 5 tine tips of one rake on the rest frame; they are kept as clicked. The view starts on the
% rake: around its tips `near` ([centerline tip] rows, 0-based) when given, else the half of the frame it
% comes from. The scroll wheel zooms where the pointer is; right-click or Backspace removes the last click;
% Done (or Enter) takes the 5 clicks. Closing the window cancels the run.
if ~usejava('desktop')
    error('runDIC:tines', 'Clicking the tine tips needs the MATLAB desktop (this MATLAB runs without one).');
end
if nargin < 5, near = []; end
[H, W] = size(gray8);
fig = figure('Name', sprintf('Click the 5 tine tips of the %s rake', upper(name)), 'NumberTitle', 'off', ...
             'Color', 'w', 'WindowState', 'maximized', 'MenuBar', 'none', 'ToolBar', 'none', 'Pointer', 'crosshair', ...
             'WindowScrollWheelFcn', @onScroll, 'KeyPressFcn', @onKey, 'CloseRequestFcn', @(~, ~) finish(false));
ax = axes(fig, 'Units', 'normalized', 'Position', [0.02 0.09 0.96 0.83]);
hImg = imagesc(ax, gray8, [0 255]); colormap(ax, gray); axis(ax, 'image', 'off'); hold(ax, 'on');
disableDefaultInteractivity(ax); ax.Toolbar.Visible = 'off';
set(hImg, 'ButtonDownFcn', @onClick);
hMarks = plot(ax, NaN, NaN, 'r+', 'MarkerSize', 16, 'LineWidth', 1.5, 'PickableParts', 'none');
head = sprintf('Click the tip of each of the 5 tines of the %s rake (the end of the rod, where it meets the tissue)', upper(name));
hs = util.buttonRow(fig, {'Whole image', 'Back to the rake', 'Remove last click', 'Done'}, ...
                    {@(~, ~) set(ax, 'XLim', [0.5 W + 0.5], 'YLim', [0.5 H + 0.5]), @(~, ~) rakeView(), ...
                     @(~, ~) undo(), @(~, ~) finish(true)}, 0, 14);   % buttons as in MATLAB's own dialogs
btnDone = hs(4); set(btnDone, 'Enable', 'off');
set(fig, 'SizeChangedFcn', @(f, ~) placeRow(f, ax, hs, 0.07));
placeRow(fig, ax, hs, 0.07);
x = zeros(0, 1); y = x; ok = false;
rakeView(); update();
uiwait(fig);
if isvalid(fig), delete(fig); end
if ~ok, util.cancel(sprintf('clicking the %s rake', name)); end
if hor, c = y - 1; t = x - 1; else, c = x - 1; t = y - 1; end     % 0-based (centerline, tip)
[c, o] = sort(c); t = t(o);
R = struct('horizontal', hor, 'out', out, 'tips', round([c(:), t(:)]));

    function rakeView()
        if ~isempty(near)
            if hor, px = near(:, 2) + 1; py = near(:, 1) + 1; else, px = near(:, 1) + 1; py = near(:, 2) + 1; end
            pad = max(60, 0.5 * max(max(px) - min(px), max(py) - min(py)));
            xs = [min(px) - pad, max(px) + pad]; ys = [min(py) - pad, max(py) + pad];
        else                                         % the half of the frame the rake comes from
            xs = [0.5 W + 0.5]; ys = [0.5 H + 0.5];
            if hor, xs = util.ternary(out < 0, [0.5 W / 2], [W / 2 W + 0.5]);
            else, ys = util.ternary(out < 0, [0.5 H / 2], [H / 2 H + 0.5]); end
        end
        set(ax, 'XLim', [max(xs(1), 0.5) min(xs(2), W + 0.5)], 'YLim', [max(ys(1), 0.5) min(ys(2), H + 0.5)]);
    end

    function onClick(~, evt)
        p = evt.IntersectionPoint;
        if evt.Button == 3, undo(); return; end
        if evt.Button ~= 1 || numel(x) >= 5, return; end
        x(end + 1, 1) = p(1); y(end + 1, 1) = p(2);
        update();
    end

    function undo()
        if isempty(x), return; end
        x(end) = []; y(end) = [];
        update();
    end

    function update()
        set(hMarks, 'XData', x, 'YData', y);
        set(btnDone, 'Enable', util.ternary(numel(x) == 5, 'on', 'off'));
        if numel(x) == 5, step = '5 of 5 clicked: Done (or Enter) takes them'; else, step = sprintf('%d of 5 clicked', numel(x)); end
        title(ax, {head, sprintf('%s. Scroll to zoom; right-click or Backspace removes the last click', step)}, 'FontSize', 11);
    end

    function onKey(~, evt)
        switch evt.Key
            case {'backspace', 'delete'}, undo();
            case {'return', 'enter'}, if numel(x) == 5, finish(true); end
        end
    end

    function onScroll(~, evt)
        p = ax.CurrentPoint(1, 1:2); xl = xlim(ax); yl = ylim(ax);
        if p(1) < xl(1) || p(1) > xl(2) || p(2) < yl(1) || p(2) > yl(2), return; end
        f = 1.25 ^ evt.VerticalScrollCount;
        wNew = min(W, diff(xl) * f); hNew = min(H, diff(yl) * f);
        if wNew < 40, return; end
        nx = p(1) - (p(1) - xl(1)) / diff(xl) * wNew; ny = p(2) - (p(2) - yl(1)) / diff(yl) * hNew;
        nx = min(max(nx, 0.5), W + 0.5 - wNew); ny = min(max(ny, 0.5), H + 0.5 - hNew);
        set(ax, 'XLim', [nx nx + wNew], 'YLim', [ny ny + hNew]);
    end

    function finish(taken)
        ok = taken;
        uiresume(fig);
    end
end


function answer = askTips(gray8, rakes, tipsPx, rakeNames, saved)
% The tips over the rest frame, in one window with the question and its answers: 'Continue',
% 'Click <RAKE> again' (that rake is then clicked), 'Skip' (with tips saved before: keep them) or
% 'Cancel' (also when the window is closed).
fig = overlay(gray8, rakes, tipsPx, '');
set(fig, 'UserData', 'Cancel', 'CloseRequestFcn', @(f, ~) util.answer(f, 'Cancel'));
ax = get(fig, 'CurrentAxes');
title(ax, {'Is every red circle on the tip of a tine? (The yellow rectangle joins the tips of each rake.)', ...
      ['If a rake''s circles are off, click that rake again with its button below.', ...
       repmat(' Skip keeps the tips saved before.', 1, saved)]}, 'FontSize', 11);
labels = [clickLabels(rakeNames), {'Continue'}, repmat({'Skip'}, 1, saved), {'Cancel'}];
hs = util.buttonRow(fig, labels, cellfun(@(a) @(~, ~) util.answer(fig, a), labels, 'UniformOutput', false), 0, 14);
for k = 1:numel(rakeNames), set(hs(k), 'Position', get(hs(k), 'Position') - [48 0 0 0]); end   % a gap between the rakes and the answers
set(fig, 'SizeChangedFcn', @(f, ~) placeRow(f, ax, hs, 0.08));
placeRow(fig, ax, hs, 0.08);
uiwait(fig);
answer = 'Cancel';
if isvalid(fig), answer = get(fig, 'UserData'); delete(fig); end
end


function placeRow(fig, ax, hs, top)
% The image above a row of buttons (pixels) centered under it, at the window's current size: on opening
% (a maximized window gets its size only once shown) and whenever it is resized.
fp = getpixelposition(fig);
set(ax, 'Position', [0.02 50 / fp(4) 0.96 1 - 50 / fp(4) - top]);
util.centerRow(hs, fp(3) / 2);
end


function ok = hasTips(testRoot)
% True when the test's test_config.json already has tine tips.
f = fullfile(testRoot, 'test_config.json');
ok = isfile(f) && isfield(util.readJSON(f), 'tine_tips_px');
end


function labels = clickLabels(rakeNames)
labels = strcat('Click', {' '}, upper(rakeNames), ' again');
end


function fig = overlay(gray8, rakes, tipsPx, path)
% Rest frame with a circle on every tip found and the tine-tip rectangle: saved to `path`, or, with an
% empty path, shown in a window (returned).
show = isempty(path);
if show
    fig = figure('Name', 'runDIC: tine tips', 'NumberTitle', 'off', 'Color', 'w', 'WindowState', 'maximized', ...
                 'MenuBar', 'none', 'ToolBar', 'none');
else
    [d, ~] = fileparts(path); if ~isfolder(d), mkdir(d); end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 900]);
end
imagesc(gray8, [0 255]); colormap(gray); axis image off; hold on
for s = fieldnames(rakes)'
    R = rakes.(s{1});
    for k = 1:size(R.tips, 1)
        if R.horizontal, x = R.tips(k, 2); y = R.tips(k, 1); else, x = R.tips(k, 1); y = R.tips(k, 2); end
        plot(x + 1, y + 1, 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
    end
end
if ~isempty(tipsPx)
    x = tipsPx.z + 1; y = tipsPx.theta + 1;
    plot([x(1) x(2) x(2) x(1) x(1)], [y(1) y(1) y(2) y(2) y(1)], 'y-', 'LineWidth', 1.5);
    pad = 150; xlim([x(1) - pad, x(2) + pad]); ylim([y(1) - pad, y(2) + pad]);
end
title('Tine tips (red circles) and tine-tip lines (yellow): every circle should sit on a tine tip', 'FontSize', 11);
if show
    drawnow
else
    print(fig, path, '-dpng', '-r100');
    close(fig);
end
end
