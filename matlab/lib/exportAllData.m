function exportAllData(testRoot, settings)
% The loading branches of every tracked deformation (settings.track_deformations, Def01-Def19) in the
% layout of CellScale's AllData export: stress_stretch/stress-stretch test<N>.AllData.xlsx, one row per
% image frame. As in CellScale's export, everything is relative to the first frame of each cycle (stretch
% 1, forces tared there). The first 27 columns are CellScale's: the Data.csv columns, CellScale's 4 points
% X1...Y4 (Coords1.csv), then lx txx ly tyy ShearAngleDeg kx txy ky tyx (allDataColumns) from the DIC:
% the homogeneous deformation gradient of every point tracked in both frames, the tine spans as widths
% and the thickness of util.thickness. The same columns follow with the suffix _4pt from CellScale's 4
% points, its rake separations at the start of the test as widths and the same thickness: CellScale's
% own values, for comparison.

cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
w = cfg.widths_mm.tine_span;                       % [w_theta w_z]: across the X force, across the Y force
[H, Hsource] = util.thickness(testRoot, settings);
P = util.phaseTiming(testRoot);
D = util.dataTable(testRoot);
raw = readtable(util.firstFile(testRoot, '*Data.csv'), 'VariableNamingRule', 'preserve', 'Encoding', 'latin1', ...
                'Delimiter', ',', 'TextType', 'char');
[~, ia] = unique(raw.Time_S, 'first'); raw = raw(sort(ia), :);
rawT = round(raw.Time_S * 10) / 10;
temp = raw{:, find(startsWith(raw.Properties.VariableNames, 'Temperature'), 1)};
co = readtable(util.firstFile(testRoot, 'Coords*.csv'), 'NumHeaderLines', 1, 'VariableNamingRule', 'preserve', 'Delimiter', ',');
[coT, ic] = unique(round(co.("Time(S)") * 10) / 10, 'first');
XY = [co.X1 co.Y1 co.X2 co.Y2 co.X3 co.Y3 co.X4 co.Y4]; XY = XY(ic, :);
wRake = [raw.Y_Size_um(1), raw.X_Size_um(1)] / 1000;   % CellScale's widths: rake separations at the start
dN = 0.25 * [-1 -1; 1 -1; 1 1; -1 1];
corners = @(r) [r(1) r(2); r(5) r(6); r(7) r(8); r(3) r(4)];   % points 1, 3, 4, 2: TL, TR, BR, BL

out = {}; done = []; missing = []; dropped = {};
for n = settings.track_deformations
    fd = util.defFolder(testRoot, n);
    f = fullfile(testRoot, 'DIC_results', [fd '.mat']);
    if ~isfile(f), missing(end + 1) = n; continue; end %#ok<AGROW>
    trk = load(f, 'times', 'posx', 'posy', 'valid');
    [t0, t1] = util.loadingBranch(testRoot, n, P, D);
    k = find(trk.times >= t0 - 1e-6 & trk.times <= t1 + 1e-6);
    t = round(trk.times(k) * 10) / 10;
    [okD, iD] = ismember(t, rawT); [okC, iC] = ismember(t, coT);
    if any(~okD), dropped{end + 1} = sprintf('Def%02d: %d frame(s) without a Data.csv row left out', n, sum(~okD)); end %#ok<AGROW>
    if any(okD & ~okC), dropped{end + 1} = sprintf('Def%02d: %d frame(s) without a Coords1.csv row (X1...Y4, _4pt empty)', n, sum(okD & ~okC)); end %#ok<AGROW>
    k = k(okD); iD = iD(okD); okC = okC(okD); iC = iC(okD);
    if isempty(k), continue; end
    m = numel(k); k0 = k(1);
    Fd = nan(2, 2, m); Fc = nan(2, 2, m);
    X0 = [trk.posx(k0, :)', trk.posy(k0, :)'];
    for q = 1:m
        v = trk.valid(k0, :)' & trk.valid(k(q), :)';
        Fd(:, :, q) = dicCore.affineFit(X0(v, :), [trk.posx(k(q), v)', trk.posy(k(q), v)']);
        if okC(q) && okC(1)
            Fc(:, :, q) = (corners(XY(iC(q), :))' * dN) / (corners(XY(iC(1), :))' * dN);
        end
    end
    dFx = raw.X_Force_mN(iD) - raw.X_Force_mN(iD(1)); dFy = raw.Y_Force_mN(iD) - raw.Y_Force_mN(iD(1));
    Cd = allDataColumns(Fd, dFx, dFy, w(1), w(2), H);
    Cc = allDataColumns(Fc, dFx, dFy, wRake(1), wRake(2), H);
    Cc.Properties.VariableNames = strcat(Cc.Properties.VariableNames, '_4pt');
    xy = nan(m, 8); xy(okC, :) = XY(iC(okC), :);
    Ti = table(raw.SetName(iD), raw.Cycle(iD), raw.Time_S(iD), raw.X_Size_um(iD), raw.Y_Size_um(iD), ...
               raw.X_Displacement_um(iD), raw.Y_Displacement_um(iD), raw.X_Force_mN(iD), raw.Y_Force_mN(iD), temp(iD), ...
               'VariableNames', {'SetName', 'Cycle', 'Time_S', 'XSize_um', 'YSize_um', 'XDisplacement_um', ...
               'YDisplacement_um', 'XForce_mN', 'YForce_mN', 'Temperature'});
    Ti = [Ti, array2table(xy, 'VariableNames', {'X1', 'Y1', 'X2', 'Y2', 'X3', 'Y3', 'X4', 'Y4'}), Cd, Cc]; %#ok<AGROW>
    out{end + 1} = Ti; done(end + 1) = n; %#ok<AGROW>
end
if isempty(out), util.log('  AllData: no tracked deformations, nothing written'); return; end
T = vertcat(out{:});
outDir = fullfile(testRoot, 'stress_stretch');
if ~isfolder(outDir), mkdir(outDir); end
[~, tn] = fileparts(testRoot);
file = fullfile(outDir, sprintf('stress-stretch %s.AllData.xlsx', regexp(tn, '^test\d+', 'match', 'once')));
if isfile(file), delete(file); end                 % writetable would keep rows of a longer old sheet
try
    writetable(T, file, 'Sheet', 'Sheet1');
catch err
    util.log('  AllData: could not write %s (open in Excel?): %s', file, err.message); return
end
fid = fopen(fullfile(outDir, 'README.txt'), 'a');
fprintf(fid, ['\n%s: the loading branches of Def%02d-Def%02d in the layout of CellScale''s AllData export, one row per ' ...
              'image frame, relative to the first frame of each cycle; lx ... tyx from the DIC (tine spans, H above), ' ...
              'the same with _4pt from CellScale''s 4 points (rake separations %.2f / %.2f mm, H above).\n'], ...
        file(numel(outDir) + 2:end), min(done), max(done), wRake(1), wRake(2));
fclose(fid);
msg = sprintf('  AllData: %d deformations (Def%02d-Def%02d), %d rows -> %s', numel(done), min(done), max(done), height(T), ...
              file(numel(testRoot) + 2:end));
if ~isempty(missing), msg = sprintf('%s; not tracked, left out: %s', msg, strjoin(compose('Def%02d', missing), ', ')); end
if ~isfinite(H), msg = [msg '; NO THICKNESS: stresses NaN']; else, msg = sprintf('%s (H = %.3f mm, %s)', msg, H, Hsource); end
util.log('%s', msg);
for q = 1:numel(dropped), util.log('    %s', dropped{q}); end
end
