function T = exportStressStretch(testRoot, settings)
% Stress-stretch data of the loading branches (Stretch start to the rake peak, or as set by hand in
% loading_branches.json) of the deformations settings.export_deformations (Def11-Def19 by default), relative to
% the common reference (first frame of the reference deformation):
%   lambda_z, lambda_theta   stretches from the DIC field: the homogeneous deformation gradient F fitted by
%                            least squares to every valid tracked point (lambda = sqrt of the diagonal of F'F)
%   force_*_mN               rake forces, tared at the start of each loading branch
%   tension_*_N_per_m        membrane tension (force per current width): F_z lambda_z / w_theta, F_theta lambda_theta / w_z
%   stress_*_kPa             Cauchy stress, tension / thickness H (NaN when no thickness is set)
%   lambda_*_cellscale       CellScale's own 4-point tracking (Coords1.csv), for comparison
%   lambda_*_rake            rake separation (X/Y_Size_um), for comparison
% Widths w_theta, w_z are the tine spans (measureTineSpans). Writes, as settings.export_formats asks,
% stress_stretch/<protocol>.csv ('csv'), stress_stretch/all_protocols.csv ('all_csv') and the plot
% stress_stretch/stress_stretch.png, .pdf, .fig ('png', 'pdf', 'fig'); and always stress_stretch/README.txt.
% With 'precond', the tracked preconditioning cycles (Def01-Def10) get a plot of their own,
% stress_stretch/stress_stretch_preconditioning.* (in the plot formats picked, else PNG), and are left out of
% the main plot.

cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
w = cfg.widths_mm.tine_span;                   % [w_theta w_z] mm
[H, Hsource] = util.thickness(testRoot);
P = util.phaseTiming(testRoot);
D = util.dataTable(testRoot);
C4 = cellscalePoints(testRoot);
outDir = fullfile(testRoot, 'stress_stretch');
if ~isfolder(outDir), mkdir(outDir); end

T = table();
tRef = [];
handSet = {};
defs = settings.export_deformations; formats = settings.export_formats;
pre = [];                                      % preconditioning cycles for their own plot
if any(strcmp(formats, 'precond'))
    d = dir(fullfile(testRoot, 'DIC_results', 'Def*.mat'));
    pre = intersect(1:10, cellfun(@(x) sscanf(x, 'Def%d'), {d.name}));
    if isempty(pre), util.log('  preconditioning plot: Def01-Def10 are not tracked'); end
end
if isempty(defs) && isempty(pre), util.log('  stress-stretch: no deformation to export'); return; end
for n = union(defs, pre)
    fd = util.defFolder(testRoot, n);
    f = fullfile(testRoot, 'DIC_results', [fd '.mat']);
    if ~isfile(f), error('runDIC:missing', '%s is not tracked (no DIC_results/%s.mat): run the track step', fd, fd); end
    trk = load(f);
    if isempty(tRef)
        tok = regexp(trk.reference, 't=([\d.]+)s', 'tokens', 'once');
        tRef = str2double(tok{1});
        P0 = c4At(C4, tRef);
        size0 = [util.atTimes(D, 'X_Size_um', tRef), util.atTimes(D, 'Y_Size_um', tRef)];
    end
    [t0, t1, hand] = util.loadingBranch(testRoot, n, P, D);
    if hand, handSet{end + 1} = sprintf('Def%02d %.1f-%.1f s', n, t0, t1); end %#ok<AGROW>
    if t1 > trk.times(end) + 0.05
        util.log('  %s: the loading branch ends at %.1f s, after the last tracked frame (%.1f s): run the track step again', ...
                 fd, t1, trk.times(end));
    end
    k = find(trk.times >= t0 - 1e-6 & trk.times <= t1 + 1e-6);
    t = trk.times(k);
    fz = util.atTimes(D, 'X_Force_mN', t); ft = util.atTimes(D, 'Y_Force_mN', t);
    fz = fz - fz(1); ft = ft - ft(1);
    lz = nan(numel(k), 1); lt = lz; lzc = lz; ltc = lz; vp = lz;
    for q = 1:numel(k)
        v = trk.valid(k(q), :)';
        x = [trk.posx(k(q), :)', trk.posy(k(q), :)'];
        F = dicCore.affineFit(trk.X0(v, :), x(v, :));
        Cg = F' * F; lz(q) = sqrt(Cg(1, 1)); lt(q) = sqrt(Cg(2, 2));
        vp(q) = 100 * mean(v);
        Pq = c4At(C4, t(q));
        if ~isempty(P0) && all(isfinite(Pq(:)))
            Fc = dicCore.bilinearCenterF(P0, Pq);
            Cc = Fc' * Fc; lzc(q) = sqrt(Cc(1, 1)); ltc(q) = sqrt(Cc(2, 2));
        end
    end
    lzr = util.atTimes(D, 'X_Size_um', t) / size0(1);
    ltr = util.atTimes(D, 'Y_Size_um', t) / size0(2);
    Tz = fz .* lz / w(1); Tt = ft .* lt / w(2);                % mN/mm = N/m
    name = regexprep(fd, '^Def\d+_', '');
    Ti = table(repmat({sprintf('Def%02d', n)}, numel(k), 1), repmat({name}, numel(k), 1), t, fz, ft, lz, lt, Tz, Tt, ...
               Tz / H, Tt / H, lzc, ltc, lzr, ltr, vp, ...
               'VariableNames', {'protocol', 'protocol_name', 'time_s', 'force_z_mN', 'force_theta_mN', ...
               'lambda_z', 'lambda_theta', 'tension_z_N_per_m', 'tension_theta_N_per_m', 'stress_zz_kPa', ...
               'stress_thetatheta_kPa', 'lambda_z_cellscale', 'lambda_theta_cellscale', 'lambda_z_rake', ...
               'lambda_theta_rake', 'valid_points_pct'});
    if ismember(n, defs) && any(strcmp(formats, 'csv')), writetable(Ti, fullfile(outDir, [fd '.csv'])); end
    T = [T; Ti]; %#ok<AGROW>
end
defName = @(v) arrayfun(@(n) sprintf('Def%02d', n), v, 'UniformOutput', false);
Tpre = T(ismember(T.protocol, defName(pre)), :);
T = T(ismember(T.protocol, defName(defs)), :);
if any(strcmp(formats, 'all_csv')), writetable(T, fullfile(outDir, 'all_protocols.csv')); end
writeReadme(fullfile(outDir, 'README.txt'), cfg, w, H, Hsource, trk, handSet, settings, pre);
plots = intersect({'png', 'pdf', 'fig'}, formats, 'stable');
main = T(~ismember(T.protocol, defName(pre)), :);
if ~isempty(plots) && height(main) > 0, plotStressStretch(main, H, cfg, fullfile(outDir, 'stress_stretch'), plots, ''); end
if ~isempty(pre)
    plotStressStretch(Tpre, H, cfg, fullfile(outDir, 'stress_stretch_preconditioning'), ...
                      util.ternary(isempty(plots), {'png'}, plots), 'preconditioning');
    util.log('  preconditioning plot: %s -> stress_stretch/stress_stretch_preconditioning', util.defList(pre));
end
if isfinite(H)
    util.log('  stress-stretch: %d protocols, %d rows (H = %.3f mm from %s, widths %.2f / %.2f mm) -> stress_stretch/', ...
             numel(defs), height(T), H, Hsource, w(1), w(2));
else
    util.log(['  stress-stretch: %d protocols, %d rows -> stress_stretch/ (NO THICKNESS: the stress columns are NaN; ' ...
              'run the thickness step, see README.md)'], numel(defs), height(T));
end
end


function C4 = cellscalePoints(testRoot)
% CellScale 4-point tracking: times and the 4 points (TL, TR, BR, BL) per row, from Coords1.csv.
C4 = [];
d = dir(fullfile(testRoot, 'Coords*.csv'));
if isempty(d), return; end
T = readtable(fullfile(d(1).folder, d(1).name), 'NumHeaderLines', 1, 'VariableNamingRule', 'preserve', 'Delimiter', ',');
t = round(T.("Time(S)") * 10) / 10;
[t, ia] = unique(t, 'first');
pts = [T.X1 T.Y1 T.X3 T.Y3 T.X4 T.Y4 T.X2 T.Y2];       % TL, TR, BR, BL = points 1, 3, 4, 2
C4 = struct('t', t, 'pts', pts(ia, :));
end


function P = c4At(C4, t)
% The 4 CellScale points (4 x 2) at time t, or [] / NaN when not available.
if isempty(C4), P = []; return; end
i = find(abs(C4.t - round(t * 10) / 10) < 1e-6, 1);
if isempty(i), P = nan(4, 2); else, P = reshape(C4.pts(i, :), 2, 4)'; end
end


function writeReadme(path, cfg, w, H, Hsource, trk, handSet, settings, pre)
fid = fopen(path, 'w');
fprintf(fid, 'Stress-stretch data, %s\n\n', cfg.sample);
fprintf(fid, 'One row per image frame of the loading branch of each protocol (Stretch start to the rake peak).\n');
fprintf(fid, 'Protocols exported: %s.\n', util.defList(settings.export_deformations));
if ~isempty(pre)
    fprintf(fid, 'Preconditioning cycles plotted on their own (stress_stretch_preconditioning): %s.\n', util.defList(pre));
end
TR = trackingStats(fileparts(fileparts(path)), settings, true);
fprintf(fid, '%s.\n', TR(1).verdict);
if ~isempty(handSet)
    fprintf(fid, 'Loading branches set by hand (loading_branches.json): %s.\n', strjoin(handSet, ', '));
end
fprintf(fid, 'Stretches are relative to the common reference, %s.\n\n', trk.reference);
fprintf(fid, 'lambda_z, lambda_theta      stretches from the DIC field (homogeneous F fitted to all valid points)\n');
fprintf(fid, 'force_z_mN, force_theta_mN  rake forces, tared at the start of each loading branch\n');
fprintf(fid, 'tension_*_N_per_m           membrane tension: F_z*lambda_z/w_theta, F_theta*lambda_theta/w_z\n');
fprintf(fid, 'stress_*_kPa                Cauchy stress = tension / thickness H\n');
fprintf(fid, 'lambda_*_cellscale          CellScale 4-point tracking (Coords1.csv), for comparison\n');
fprintf(fid, 'lambda_*_rake               rake separation, for comparison\n');
fprintf(fid, 'valid_points_pct            DIC points tracked in this frame\n\n');
fprintf(fid, 'z = axial (tester X), theta = circumferential (tester Y).\n');
fprintf(fid, 'w_theta = %.2f mm (divides F_z), w_z = %.2f mm (divides F_theta): tine spans.\n', w(1), w(2));
if isfinite(H), fprintf(fid, 'Thickness H = %.4f mm (%s).\n', H, Hsource);
else, fprintf(fid, 'Thickness: NOT SET, the stress columns are NaN.\n'); end
fprintf(fid, 'Image scale %.3f px/mm; DIC subset %.2f mm, grid step %.2f mm.\n', trk.px_per_mm, trk.subset_mm, trk.step_mm);
[~, notes] = util.sampleInfo(fileparts(fileparts(path)));
if ~isempty(notes), fprintf(fid, '\nNotes:\n'); fprintf(fid, '  - %s\n', notes{:}); fprintf(fid, '\n'); end
fprintf(fid, 'Written %s by runDIC (MATLAB).\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fclose(fid);
end


function plotStressStretch(T, H, cfg, base, formats, what)
% Every protocol of T on one plot: circles axial, triangles circumferential, one color per protocol; saved as
% base.png, base.pdf and/or base.fig (formats). what ('' or e.g. 'preconditioning') is added to the title.
if isfinite(H), yz = T.stress_zz_kPa; yt = T.stress_thetatheta_kPa; yl = 'Cauchy stress (kPa)';
else, yz = T.tension_z_N_per_m; yt = T.tension_theta_N_per_m; yl = 'membrane tension (N/m)'; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 600]);   % 3:2 (landscape)
hold on
cols = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; 227 119 194; 127 127 127; 188 189 34; ...
        23 190 207] / 255;
key = {'MarkerSize', 6, 'LineWidth', 1, 'LineStyle', 'none', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w'};
plot(NaN, NaN, 'o', key{:}, 'DisplayName', 'Longitudinal (z)');                  % legend keys for the symbols
plot(NaN, NaN, '^', key{:}, 'DisplayName', 'Circumferential (\theta)');
protos = unique(T.protocol, 'stable');
names = cell(size(protos));                     % legend names: Equibiaxial 1, 1-0.75 ratio, Cycle 1, ...
for i = 1:numel(protos)
    names{i} = strrep(T.protocol_name{find(strcmp(T.protocol, protos{i}), 1)}, '_', ' ');
    cyc = regexpi(names{i}, '^Preconditioning c0*(\d+)$', 'tokens', 'once');
    if ~isempty(cyc), names{i} = ['Cycle ' cyc{1}]; end
end
for u = unique(names)'                           % a name used more than once is numbered in order
    k = find(strcmp(names, u{1}));
    if numel(k) > 1, for j = 1:numel(k), names{k(j)} = sprintf('%s %d', u{1}, j); end, end
end
for i = 1:numel(protos)
    m = strcmp(T.protocol, protos{i});
    c = cols(mod(i - 1, size(cols, 1)) + 1, :);
    name = names{i};
    plot(T.lambda_z(m), yz(m), 'o', 'Color', c, 'MarkerSize', 4.5, 'LineWidth', 1, 'DisplayName', name);
    plot(T.lambda_theta(m), yt(m), '^', 'Color', c, 'MarkerSize', 4.5, 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel('Stretch'); ylabel(yl); grid on; legend('Location', 'northwest');
xlim([1 1.6]);
set(gca, 'FontSize', get(gca, 'FontSize') + 1);                 % labels, ticks, legend and title one point larger
if isfinite(H), ylim([0 120]); end                          % kPa; tensions (no thickness) keep their own range
ttl = regexprep(strrep(cfg.sample, '_', ' '), {',\s*labeled\s+', '\s*\(\d+/\d+/\d+\)'}, {', ', ''}, 'ignorecase');
if ~isempty(what), ttl = sprintf('%s: %s', ttl, what); end
title(ttl, 'Interpreter', 'none');
for f = formats
    switch f{1}
        case 'png', print(fig, [base '.png'], '-dpng', '-r100');
        case 'pdf', exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
        case 'fig'
            set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');   % drawn hidden; opens shown
            savefig(fig, [base '.fig']);
    end
end
close(fig);
end
