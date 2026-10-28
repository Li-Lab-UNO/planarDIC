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
tRec = arrayfun(@(i) recoverStart(P, i), 1:numel(names));
checkSort(testRoot, imgDir, names, t0, t1, tRec, nUnassigned, settings);
end


function t = recoverStart(P, n)
% Where the n-th cycle's loading ends and unloading (its Recover phase) starts (s).
[~, t] = util.stretchPhase(P, n);
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


function checkSort(testRoot, imgDir, names, t0, t1, tRec, nUnassigned, settings)
% Completeness of every cycle, from the tester's own files: images, Data.csv rows and Coords1.csv rows
% missing (gaps in their times), and images without a Data.csv or Coords1.csv row at their time (a frame
% needs all three). The note says where the missing frames are and whether they fall in loading or
% unloading (from tRec, the Recover start). Asks whether to proceed (desktop only).
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
    miss = gaps(t, t0(i), tEnd, dt);
    missD = gaps(tD(tD >= t0(i) - 1e-6 & tD < tEnd - 1e-6), t0(i), tEnd, dtD);
    missing(i) = numel(miss); dataMissing(i) = numel(missD);
    missC = [];
    if ~isempty(tC)
        missC = gaps(tC(tC >= t0(i) - 1e-6 & tC < tEnd - 1e-6), t0(i), tEnd, dtC);
        coordsMissing(i) = numel(missC);
    end
    tr = round(t * 10) / 10;
    unmatched(i) = sum(~ismember(tr, tD) | (~isempty(tC) & ~ismember(tr, tC)));
    same = @(m) numel(m) == numel(miss) && all(abs(m - miss) < 0.05);
    p = {};
    if missing(i) > 0, p{end + 1} = sprintf('%d image(s) missing%s', missing(i), where(miss, dt, tRec(i))); end %#ok<AGROW>
    rows = {missD, 'Data.csv', dtD; missC, 'Coords1.csv', dtC};
    for k = 1:2
        m = rows{k, 1};
        if isempty(m), continue; end
        if missing(i) > 0 && same(m), at = ' (the same frames)'; else, at = where(m, rows{k, 3}, tRec(i)); end
        p{end + 1} = sprintf('%d %s row(s) missing%s', numel(m), rows{k, 2}, at); %#ok<AGROW>
    end
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
             'ToolBar', 'none', 'Visible', 'off', 'Resize', 'off', ...
             'CloseRequestFcn', @(f, ~) util.answer(f, 'Cancel the run'));
rows = [T.folder, num2cell(T.images), compose('%.1f', T.first_image_s), compose('%.1f', T.last_image_s), ...
        num2cell([T.images_missing, T.data_rows_missing, T.coords_rows_missing]), T.note];
rows(flag, end) = strcat('<html><b><font color="red">', rows(flag, end), '</font></b></html>');
cw = {190, 60, 75, 75, 70, 100, 100, 340};
tbl = uitable(fig, 'Data', rows, 'ColumnName', {'Folder', 'Images', 'First|image (s)', 'Last|image (s)', 'Images|missing', ...
        'Data.csv rows|missing', 'Coords1.csv|rows missing', 'Note'}, 'ColumnWidth', cw, 'RowName', []);
tw = sum([cw{:}]) + 3; th = 44 + 23 * n;             % the table's own size (two-line header, rows): no empty space around it
W = tw + 40;
[~, tn] = fileparts(testRoot);
info = [{sprintf(['%s: %d images, one every %.2f s; %d images outside every cycle (before the first or after ' ...
        'the last).'], strrep(tn, '_', ' '), sum(count), dt, nUnassigned)}; strcat(names(flag)', {': '}, note(flag))];
if any(flag)
    q = ['Some folders are missing frames (red). A few missing frames are fine (the camera drops frames now and then), ' ...
         'and frames missing during unloading don''t affect the loading curves; a folder with no images or many missing ' ...
         'means the test or the folders need a look. Proceed with this test?'];
else
    q = 'Every folder has all its images, Data.csv rows and Coords1.csv rows. Proceed with this test?';
end
hInfo = uicontrol(fig, 'Style', 'text', 'Position', [20 0 W - 40 20], 'BackgroundColor', 'w', ...
                  'HorizontalAlignment', 'left', 'FontSize', 10);   % the notes in full, under the table
hQ = uicontrol(fig, 'Style', 'text', 'Position', [20 0 W - 40 20], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
               'FontSize', 10, 'FontWeight', 'bold');
hi = fitText(hInfo, info); hq = fitText(hQ, {q});
yq = 14 + 22 + 12; yi = yq + hq + 10; yt = yi + hi + 12; H = yt + th + 16;   % from the bottom: buttons, question, notes, table
set(hQ, 'Position', [20 yq W - 40 hq]); set(hInfo, 'Position', [20 yi W - 40 hi]);
set(tbl, 'Position', [20 yt tw th]);
set(fig, 'Position', [100 100 W H]); movegui(fig, 'center');
answers = {'Yes, proceed', 'Cancel the run'};
if isfield(settings, 'n_tests') && settings.n_tests > 1       % other tests in the run: stopping this one is a choice
    answers = {'Yes, proceed', 'No, skip to the next test', 'Cancel the run'};
end
util.buttonRow(fig, answers, cellfun(@(a) @(~, ~) util.answer(fig, a), answers, 'UniformOutput', false), W - 20, 14);
set(fig, 'Visible', 'on');
uiwait(fig);
answer = 'Cancel the run';
if isvalid(fig), answer = get(fig, 'UserData'); delete(fig); end
if strcmp(answer, 'Cancel the run'), util.cancel('the sorting check'); end
if strcmp(answer, 'No, skip to the next test')
    error('runDIC:stopped', 'Stopped after sorting (see Images/sort_summary.csv)');
end
end



function h = fitText(hText, lines)
% Wrap `lines` to the width of the text control and return the height (pixels) they take.
w = textwrap(hText, lines);
set(hText, 'String', w);
e = get(hText, 'Extent'); h = ceil(e(4)) + 2;
end


function n = cleanName(name)
% Folder-safe protocol name: 'Equibiaxial save cycle' -> 'Equibiaxial', '1:0.75 ratio' -> '1-0.75_ratio' style.
n = regexprep(name, '\s+save\s+cyc?le$', '', 'ignorecase');     % 'save cycle' (also 'save cyle'): a label only
n = regexprep(strtrim(n), '[^A-Za-z0-9.\-]+', '_');
n = regexprep(n, '^_+|_+$', '');
end


function miss = gaps(t, tStart, tEnd, dt)
% Times (s) of the frames missing from times t between tStart and tEnd, at one every dt.
u = [tStart - dt; t(:); tEnd];
k = max(0, round(diff(u) / dt) - 1);
miss = zeros(0, 1);
for i = find(k)'
    miss = [miss; u(i) + dt * (1:k(i))']; %#ok<AGROW>
end
end


function s = where(miss, dt, tRec)
% Where the missing frames are: their times, as runs, and whether they fall in loading or unloading.
runs = {}; first = 1;
for i = 2:numel(miss) + 1
    if i > numel(miss) || miss(i) - miss(i - 1) > 1.5 * dt
        if i - 1 == first, runs{end + 1} = sprintf('%.1f s', miss(first)); %#ok<AGROW>
        else, runs{end + 1} = sprintf('%.1f-%.1f s', miss(first), miss(i - 1)); end %#ok<AGROW>
        first = i;
    end
end
if numel(runs) > 3, runs = [runs(1:3), {'...'}]; end
n = sum(miss < tRec - 1e-6);
if n == numel(miss), phase = 'all loading';
elseif n == 0, phase = 'all unloading';
else, phase = sprintf('%d loading, %d unloading', n, numel(miss) - n); end
s = sprintf(': %s (%s)', strjoin(runs, ', '), phase);
end
