function sortImages(testRoot, settings)
% Sort the CellScale images of a test into one folder per load-unload cycle, Images/DefNN_<name>.
% Cycle boundaries come from Logs/*_PhaseTiming.csv: each cycle runs from its Stretch start to the next
% Stretch start (the last one to the End time). Images already in Images/ subfolders are found and moved
% too, so running it again is safe. Every move is recorded in Images/sort_manifest.csv. Afterwards every
% cycle is checked for missing images, Data.csv rows and Coords1.csv rows (checkSort); the table goes to
% Images/sort_summary.csv and the Command Window, and with settings.review_sort it is shown in a window
% that asks whether to proceed.

P = util.phaseTiming(testRoot);
stretch = P(strcmp({P.Phase}, 'Stretch'));
tEnd = P(strcmp({P.Phase}, 'End')).Start;
sets = [stretch.Set];
names = cell(1, numel(stretch)); t0 = zeros(1, numel(stretch)); t1 = t0;
for i = 1:numel(stretch)
    s = stretch(i);
    t0(i) = s.Start;
    if i < numel(stretch), t1(i) = stretch(i + 1).Start; else, t1(i) = tEnd + 1; end
    names{i} = sprintf('Def%02d_%s', i, cleanName(s.SetName));
    if sum(sets == s.Set) > 1
        names{i} = sprintf('%s_c%02d', names{i}, s.Cycle);
    end
end

imgDir = fullfile(testRoot, 'Images');
% a test sorted before (by this code or the Python pipeline) keeps its folder names, so its tracks, which
% are named after them, still match (e.g. Def19_Equibiaxial_Save from older versions)
for i = 1:numel(names)
    old = dir(fullfile(imgDir, sprintf('Def%02d_*', i)));
    old = old([old.isdir]);
    if isscalar(old), names{i} = old.name; end
end
files = dir(fullfile(imgDir, '**', '*.jpg'));
moves = {}; nUnassigned = 0;
for k = 1:numel(files)
    t = util.frameTime(files(k).name);
    if isnan(t), continue; end
    i = find(t0 - 1e-6 <= t & t < t1 - 1e-6, 1);
    if isempty(i), nUnassigned = nUnassigned + 1; continue; end
    old = fullfile(files(k).folder, files(k).name);
    new = fullfile(imgDir, names{i}, files(k).name);
    if ~strcmp(old, new), moves(end + 1, :) = {old, new}; end %#ok<AGROW>
end
util.log('Sorting %d images into %d deformation folders (%d to move, %d outside every cycle)', ...
         numel(files), numel(names), size(moves, 1), nUnassigned);
if ~isempty(moves), moveImages(moves, testRoot, imgDir); end
checkSort(testRoot, imgDir, names, t0, t1, nUnassigned, settings);
end


function moveImages(moves, testRoot, imgDir)
man = fullfile(imgDir, 'sort_manifest.csv');
newManifest = ~isfile(man);
fid = fopen(man, 'a');
if newManifest, fprintf(fid, 'old_path,new_path\n'); end
rel = @(p) strrep(p(numel(testRoot) + 2:end), '\', '/');
oldDirs = {};
for k = 1:size(moves, 1)
    [d, ~] = fileparts(moves{k, 2});
    if ~isfolder(d), mkdir(d); end
    if isfile(moves{k, 2}), fclose(fid); error('runDIC:exists', 'Already exists: %s', moves{k, 2}); end
    movefile(moves{k, 1}, moves{k, 2});
    fprintf(fid, '"%s","%s"\n', rel(moves{k, 1}), rel(moves{k, 2}));
    oldDirs{end + 1} = fileparts(moves{k, 1}); %#ok<AGROW>
end
fclose(fid);
% remove folders the images were moved out of, if now empty
for d = unique(oldDirs)
    if ~strcmp(d{1}, imgDir) && isfolder(d{1})
        left = dir(d{1}); left = left(~ismember({left.name}, {'.', '..', '.DS_Store', 'Thumbs.db'}));
        if isempty(left), rmdir(d{1}, 's'); end
    end
end
end


function checkSort(testRoot, imgDir, names, t0, t1, nUnassigned, settings)
% Completeness of every cycle, from the tester's own files: images, Data.csv rows and Coords1.csv rows
% missing (gaps in their times), and images without a Data.csv or Coords1.csv row at their time (a frame
% needs all three). Asks whether to proceed (desktop only).
n = numel(names);
times = cell(n, 1);
for i = 1:n, times{i} = util.findFrames(fullfile(imgDir, names{i})); end
dt = median(cell2mat(cellfun(@(t) diff(t(:)), times, 'UniformOutput', false)));   % image interval (s)
D = util.dataTable(testRoot); tD = D.t; dtD = median(diff(tD));
tC = []; dtC = NaN;
c = dir(fullfile(testRoot, 'Coords*.csv'));
if ~isempty(c)
    C = readtable(fullfile(c(1).folder, c(1).name), 'NumHeaderLines', 1, 'VariableNamingRule', 'preserve', 'Delimiter', ',');
    tC = unique(round(C.("Time(S)") * 10) / 10); dtC = median(diff(tC));
end
count = zeros(n, 1); first = nan(n, 1); last = first; missing = zeros(n, 1); dataMissing = missing;
coordsMissing = nan(n, 1); unmatched = missing; note = cell(n, 1);
for i = 1:n
    t = times{i}(:);
    count(i) = numel(t);
    note{i} = '';
    if isempty(t), note{i} = 'NO IMAGES'; continue; end
    first(i) = t(1); last(i) = t(end);
    tEnd = min(t1(i), t(end) + dt);
    [missing(i), gap] = gaps(t, t0(i), tEnd, dt);
    dataMissing(i) = gaps(tD(tD >= t0(i) - 1e-6 & tD < tEnd - 1e-6), t0(i), tEnd, dtD);
    if ~isempty(tC), coordsMissing(i) = gaps(tC(tC >= t0(i) - 1e-6 & tC < tEnd - 1e-6), t0(i), tEnd, dtC); end
    tr = round(t * 10) / 10;
    unmatched(i) = sum(~ismember(tr, tD) | (~isempty(tC) & ~ismember(tr, tC)));
    p = {};
    if missing(i) > 0, p{end + 1} = sprintf('%d image(s) missing (longest gap %.1f s)', missing(i), gap); end %#ok<AGROW>
    if dataMissing(i) > 0, p{end + 1} = sprintf('%d Data.csv row(s) missing', dataMissing(i)); end %#ok<AGROW>
    if coordsMissing(i) > 0, p{end + 1} = sprintf('%d Coords1.csv row(s) missing', coordsMissing(i)); end %#ok<AGROW>
    if unmatched(i) > 0, p{end + 1} = sprintf('%d image(s) without a Data/Coords row', unmatched(i)); end %#ok<AGROW>
    note{i} = strjoin(p, '; ');
end
T = table(names(:), count, round(first, 1), round(last, 1), round(t1(:) - t0(:), 1), missing, dataMissing, ...
          coordsMissing, unmatched, note, 'VariableNames', {'folder', 'images', 'first_image_s', 'last_image_s', ...
          'cycle_length_s', 'images_missing', 'data_rows_missing', 'coords_rows_missing', 'images_unmatched', 'note'});
T.cycle_length_s(end) = NaN;                                   % the last cycle runs to the End time
writetable(T, fullfile(imgDir, 'sort_summary.csv'));
flag = ~cellfun(@isempty, note);
util.log('  %d folders, %d images (one every %.2f s); %d outside every cycle; %s -> Images/sort_summary.csv', n, ...
         sum(count), dt, nUnassigned, util.ternary(any(flag), 'gaps below', 'images, Data.csv and Coords1.csv complete'));
for i = find(flag)', util.log('  %s: %s', names{i}, note{i}); end
if ~(settings.review_sort && usejava('desktop')), return; end

fig = figure('Name', 'runDIC: sorted images', 'NumberTitle', 'off', 'Color', 'w', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Position', [100 100 1040 620]);
movegui(fig, 'center');
rows = [T.folder, num2cell(T.images), compose('%.1f', T.first_image_s), compose('%.1f', T.last_image_s), ...
        num2cell([T.images_missing, T.data_rows_missing, T.coords_rows_missing]), T.note];
rows(flag, end) = strcat('<html><b><font color="red">', rows(flag, end), '</font></b></html>');
uitable(fig, 'Data', rows, 'ColumnName', {'Folder', 'Images', 'First|image (s)', 'Last|image (s)', 'Images|missing', ...
        'Data.csv rows|missing', 'Coords1.csv|rows missing', 'Note'}, 'ColumnWidth', {190, 60, 75, 75, 70, 100, 100, 280}, ...
        'RowName', [], 'Units', 'normalized', 'Position', [0.02 0.12 0.96 0.86]);
[~, tn] = fileparts(testRoot);
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.02 0.02 0.96 0.07], 'BackgroundColor', 'w', ...
          'HorizontalAlignment', 'left', 'FontSize', 11, 'String', sprintf(['%s: %d images, one every %.2f s; %d images ' ...
          'outside every cycle (before the first or after the last).'], strrep(tn, '_', ' '), sum(count), dt, nUnassigned));
if any(flag)
    q = {'Some folders are missing frames (red).', 'A few missing frames are fine (the camera drops frames now and then);', ...
         'a folder with no images or many missing means the test or the folders need a look.', '', 'Proceed with this test?'};
else
    q = {'Every folder has all its images, Data.csv rows and Coords1.csv rows.', '', 'Proceed with this test?'};
end
answer = questdlg(q, 'runDIC: sorted images', 'Yes, proceed', 'No, stop this test', 'Yes, proceed');
if isvalid(fig), close(fig); end
if strcmp(answer, 'No, stop this test')
    error('runDIC:stopped', 'Stopped after sorting (see Images/sort_summary.csv)');
end
end


function n = cleanName(name)
% Folder-safe protocol name: 'Equibiaxial save cycle' -> 'Equibiaxial', '1:0.75 ratio' -> '1-0.75_ratio' style.
n = regexprep(name, '\s+save\s+cyc?le$', '', 'ignorecase');     % 'save cycle' (also 'save cyle'): a label only
n = regexprep(strtrim(n), '[^A-Za-z0-9.\-]+', '_');
n = regexprep(n, '^_+|_+$', '');
end


function [m, longest] = gaps(t, tStart, tEnd, dt)
% Frames missing from times t (s) between tStart and tEnd at one every dt, and the longest gap.
d = diff([tStart; t(:); tEnd]);
m = sum(max(0, round(d / dt) - 1)); longest = max(d);
end
