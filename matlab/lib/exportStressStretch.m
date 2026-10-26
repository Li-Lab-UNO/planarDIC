function T = exportStressStretch(testRoot, settings)
% Stress-stretch data of the loading branches (Stretch start to the rake peak, or as set by hand in
% loading_branches.json) of the analyzed deformations (settings.deformations, Def11-Def19), relative to
% the common reference (first frame of the reference deformation):
%   lambda_z, lambda_theta   stretches from the DIC field: the homogeneous deformation gradient F fitted by
%                            least squares to every valid tracked point (lambda = sqrt of the diagonal of F'F)
%   force_*_mN               rake forces, tared at the start of each loading branch
%   tension_*_N_per_m        membrane tension (force per current width): F_z lambda_z / w_theta, F_theta lambda_theta / w_z
%   stress_*_kPa             Cauchy stress, tension / thickness H (NaN when no thickness is set)
%   lambda_*_cellscale       CellScale's own 4-point tracking (Coords1.csv), for comparison
%   lambda_*_rake            rake separation (X/Y_Size_um), for comparison
% Widths w_theta, w_z are the tine spans (measureTineSpans). Writes stress_stretch/<protocol>.csv,
% stress_stretch/all_protocols.csv, stress_stretch/README.txt and stress_stretch/stress_stretch.png.

cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
w = cfg.widths_mm.tine_span;                   % [w_theta w_z] mm
[H, Hsource] = util.thickness(testRoot, settings);
P = util.phaseTiming(testRoot);
D = util.dataTable(testRoot);
C4 = cellscalePoints(testRoot);
outDir = fullfile(testRoot, 'stress_stretch');
if ~isfolder(outDir), mkdir(outDir); end

T = table();
tRef = [];
handSet = {};
for n = settings.deformations
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
    writetable(Ti, fullfile(outDir, [fd '.csv']));
    T = [T; Ti]; %#ok<AGROW>
end
writetable(T, fullfile(outDir, 'all_protocols.csv'));
writeReadme(fullfile(outDir, 'README.txt'), cfg, w, H, Hsource, trk, handSet, settings);
plotStressStretch(T, H, cfg, fullfile(outDir, 'stress_stretch.png'));
if isfinite(H)
    util.log('Stress-stretch: %d protocols, %d rows (H = %.3f mm from %s, widths %.2f / %.2f mm) -> stress_stretch/', ...
             numel(settings.deformations), height(T), H, Hsource, w(1), w(2));
else
    util.log(['Stress-stretch: %d protocols, %d rows -> stress_stretch/ (NO THICKNESS: the stress columns are NaN; ' ...
              'run the thickness step, see README.md)'], numel(settings.deformations), height(T));
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


function writeReadme(path, cfg, w, H, Hsource, trk, handSet, settings)
fid = fopen(path, 'w');
fprintf(fid, 'Stress-stretch data, %s\n\n', cfg.sample);
fprintf(fid, 'One row per image frame of the loading branch of each protocol (Stretch start to the rake peak).\n');
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


function plotStressStretch(T, H, cfg, path)
% Every protocol on one plot: circles axial, triangles circumferential, one color per protocol.
if isfinite(H), yz = T.stress_zz_kPa; yt = T.stress_thetatheta_kPa; yl = 'Cauchy stress (kPa)';
else, yz = T.tension_z_N_per_m; yt = T.tension_theta_N_per_m; yl = 'membrane tension (N/m)'; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
hold on
cols = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; 227 119 194; 127 127 127; 188 189 34] / 255;
protos = unique(T.protocol, 'stable');
for i = 1:numel(protos)
    m = strcmp(T.protocol, protos{i});
    c = cols(mod(i - 1, 9) + 1, :);
    plot(T.lambda_z(m), yz(m), 'o', 'Color', c, 'MarkerSize', 4.5, 'LineWidth', 1, 'DisplayName', sprintf('%s %s', protos{i}, T.protocol_name{find(m, 1)}));
    plot(T.lambda_theta(m), yt(m), '^', 'Color', c, 'MarkerSize', 4.5, 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel('stretch'); ylabel(yl); grid on; legend('Location', 'northwest', 'Interpreter', 'none');
title({strrep(cfg.sample, '_', ' '), 'circles: axial (\lambda_z), triangles: circumferential (\lambda_\theta)'});
print(fig, path, '-dpng', '-r100');
close(fig);
end
