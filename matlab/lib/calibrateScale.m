function result = calibrateScale(testRoot, settings)
% Image scale (px per mm) from the rake motion. The 20 tine tips are tracked in every frame of the
% calibration deformations (Def01 and Def19 by default); the change in tip-to-tip separation (px) of each
% opposing tine pair is regressed against the change in rake separation the tester reports (X/Y_Size_um).
% Tine pairs that disagree with the rest are left out: more than 6% from the median of all pairs, the two
% deformations more than 5% apart, or not tracked throughout (calibration_exclude in test_config.json,
% e.g. ["Y3"], replaces this rule). The result is the mean of the middle three pairs of both axes.
% Writes scale_calibration.json and DIC_results/scale_calibration.png.

cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
tines = cfg.calibration_tines_px;
D = util.dataTable(testRoot);
rows = [];                  % [deformation axis(1=X,2=Y) pair px_per_mm rms_px]
plt = {};
for dn = settings.calibration_deformations
    folder = util.defFolder(testRoot, dn);
    [times, paths] = util.findFrames(fullfile(testRoot, 'Images', folder));
    util.log('  %s: tracking 20 tine tips in %d frames', folder, numel(paths));
    tips = trackTips(paths, tines);
    for ax = 1:2
        if ax == 1, lo = 'left'; hi = 'right'; size_ = 'X_Size_um'; force_ = 'X_Force_mN';
        else, lo = 'top'; hi = 'bottom'; size_ = 'Y_Size_um'; force_ = 'Y_Force_mN'; end
        sep = tips.(hi) - tips.(lo);                       % (T, 5)
        rk = util.atTimes(D, size_, times); rake = rk - rk(1);
        dsep = sep - sep(1, :);
        for pair = 1:size(sep, 2)
            y = dsep(:, pair); m = isfinite(y) & isfinite(rake);
            if sum(m) < 10, continue; end
            pf = polyfit(rake(m), y(m), 1);
            res = y(m) - polyval(pf, rake(m));
            rows(end + 1, :) = [dn, ax, pair, 1000 * pf(1), sqrt(mean(res .^ 2))]; %#ok<AGROW>
        end
        plt(end + 1, :) = {dn, ax, rake, dsep, util.atTimes(D, force_, times)}; %#ok<AGROW>
    end
end

% outlier rule (or the config's calibration_exclude)
axisName = 'XY';
names = arrayfun(@(a, p) sprintf('%s%d', axisName(a), p), rows(:, 2), rows(:, 3), 'UniformOutput', false);
excluded = {}; excludedBy = 'outlier rule';
if isfield(cfg, 'calibration_exclude')
    excluded = cellstr(cfg.calibration_exclude); excludedBy = 'calibration_exclude in test_config.json';
else
    ref = median(rows(:, 4), 'omitnan');
    nDef = numel(settings.calibration_deformations);
    for e = unique(names)'
        v = rows(strcmp(names, e{1}), 4);
        if numel(v) < nDef || any(~isfinite(v))
            why = 'not tracked in every deformation';
        elseif abs(mean(v) / ref - 1) > 0.06
            why = sprintf('%.2f px/mm, %+.1f%% from the median %.2f', mean(v), 100 * (mean(v) / ref - 1), ref);
        elseif numel(v) > 1 && (max(v) - min(v)) / ref > 0.05
            why = sprintf('deformations disagree: %s px/mm', strjoin(compose('%.2f', v), ', '));
        else
            continue
        end
        excluded{end + 1} = e{1}; %#ok<AGROW>
        util.log('  left out tine pair %s (%s)', e{1}, why);
    end
end
keep = ~ismember(names, excluded);
mid = keep & ismember(rows(:, 3), [2 3 4]);
if numel(unique(names(mid))) < 2
    mid = keep;
    util.log('  fewer than two middle pairs left: using every pair left');
end
if ~any(mid), error('runDIC:scale', 'Every tine pair is excluded: check DIC_results/config_tine_tips.png'); end
s = rows(mid, 4);
ppm = mean(s);
result = struct('px_per_mm', ppm, 'mm_per_px', 1 / ppm, 'um_per_px', 1000 / ppm, 'std_px_per_mm', std(s), ...
    'px_per_mm_X', mean(rows(mid & rows(:, 2) == 1, 4)), 'px_per_mm_Y', mean(rows(mid & rows(:, 2) == 2, 4)), ...
    'method', 'tine-tip separation vs reported rake separation (middle 3 tine pairs, all frames)', ...
    'excluded', {sort(excluded)}, 'excluded_by', excludedBy, ...
    'deformations', {arrayfun(@(d) util.defFolder(testRoot, d), settings.calibration_deformations, 'UniformOutput', false)});
util.writeJSON(fullfile(testRoot, 'scale_calibration.json'), result);
util.log('  scale: %.3f px/mm (X %.3f, Y %.3f; SD across pairs %.2f%%); left out: %s', ppm, result.px_per_mm_X, ...
         result.px_per_mm_Y, 100 * std(s) / ppm, strjoin(sort(excluded), ', '));

% diagnostic figure
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 450]);
cols = lines(5);
subplot(1, 3, 1); hold on
for q = 1:size(plt, 1)
    ls = '-'; if plt{q, 2} == 2, ls = '--'; end
    for pair = 1:5, plot(plt{q, 3} / 1000, plt{q, 4}(:, pair), ls, 'Color', cols(pair, :)); end
end
xlabel('reported change in rake separation (mm)'); ylabel('change in tine-tip separation (px)');
title('Tine tips vs rake motion (solid X, dashed Y)'); grid on
subplot(1, 3, 2); hold on
for q = 1:size(plt, 1)
    ls = '-'; if plt{q, 2} == 2, ls = '--'; end
    for pair = 1:5
        r = plt{q, 4}(:, pair) - ppm * plt{q, 3} / 1000;
        plot(plt{q, 5}, r - median(r(1:min(5, end)), 'omitnan'), ls, 'Color', cols(pair, :));
    end
end
xlabel('rake force (mN)'); ylabel('residual from one scale (px)'); title('Residual vs force'); grid on
subplot(1, 3, 3); hold on
for ax = 1:2
    q = rows(:, 2) == ax;
    scatter(rows(q, 3) + (ax - 1.5) * 0.15, rows(q, 4), 30, 'filled', 'DisplayName', sprintf('%s axis', axisName(ax)));
end
yline(ppm, 'k-', 'HandleVisibility', 'off');
xlabel('tine pair (1 = top/left ... 5 = bottom/right)'); ylabel('px per mm'); legend('Location', 'best');
title(sprintf('Scale per tine pair; line = result %.2f px/mm', ppm)); grid on
out = fullfile(testRoot, 'DIC_results', 'scale_calibration.png');
if ~isfolder(fileparts(out)), mkdir(fileparts(out)); end
print(fig, out, '-dpng', '-r90');
close(fig);
end


function tips = trackTips(paths, tines)
% Tip coordinate of every tine in every frame: struct of (T, 5) arrays per rake.
sides = fieldnames(tines);
state = struct(); tips = struct();
for k = 1:numel(sides)
    state.(sides{k}) = tines.(sides{k}).tips;
    tips.(sides{k}) = nan(numel(paths), size(tines.(sides{k}).tips, 1));
end
for f = 1:numel(paths)
    g = imageOps.gaussBlur(imageOps.readGray(paths{f}), 0.8);
    for k = 1:numel(sides)
        R = tines.(sides{k});
        for i = 1:size(state.(sides{k}), 1)
            c = state.(sides{k})(i, 1); t = state.(sides{k})(i, 2);
            for it = 1:2
                [c, tn, contrast] = calFindTip(g, c, t, R.horizontal, R.out);
                if isfinite(tn), t = tn; end
            end
            state.(sides{k})(i, :) = [c, t];
            if contrast > 30, tips.(sides{k})(f, i) = t; end
        end
    end
end
end
