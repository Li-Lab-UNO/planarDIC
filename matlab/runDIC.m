function runDIC(folder, varargin)
% runDIC  Stress-stretch data from a CellScale BioTester planar biaxial test.
%
%   runDIC              asks for a folder (a test folder test1_..., or a sample folder holding tests),
%                       then which steps to run
%   runDIC(folder)      runs that test, or every test in that sample folder
%   runDIC(folder, 'steps', {'track', 'export'})   only some steps
%
% Settings are in DIC_settings.m. A test folder holds the CellScale output: *Data.csv, Coords1.csv,
% Logs/*_PhaseTiming.csv and the images in Images/. The steps:
%   info    asks for the sample label (e.g. healthy) and notes -> test_config.json
%   thickness  opens the thickness window: each test's thickness typed (e.g. from the workbook) or
%           measured from its photos, filled in with the saved one -> test_config.json
%   sort    images into one folder per load-unload cycle (Images/Def01_... to Def19_...), then checks that
%           every cycle has all its images, Data.csv rows and Coords1.csv rows -> Images/sort_summary.csv
%   tines   finds the 20 tine tips in the first image of Def11 -> test_config.json
%           check DIC_results/config_tine_tips.png: every red circle should sit on a tine tip
%   scale   image scale from the tine motion in Def01 and Def19 -> scale_calibration.json
%           check DIC_results/scale_calibration.png
%   widths  tine spans (specimen widths for the stresses) -> test_config.json
%   track   DIC tracking of Def01-Def19 (the loading branches) -> DIC_results/DefNN_*.mat
%   artifacts  finds local tracking errors (a speckle particle sliding over the tissue, a shifting glare
%           patch), excludes them in test_config.json and tracks again -> DIC_results/artifact_check.png
%   export  stress-stretch CSVs -> stress_stretch/ (one per protocol and all_protocols.csv), and the loading
%           branches of Def01-Def19 in CellScale's AllData layout -> stress_stretch/stress-stretch test<N>.AllData.xlsx
%   report  the tracking report (viewTracking) of the existing tracks; it also opens at the end of a run
%           that tracked

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
s = DIC_settings();
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k + 1}; end
interactive = nargin < 1 || isempty(folder);
if interactive
    start = pwd;                     % open next to the folder of the last run, so it and the tests
    last = getpref('runDIC', 'last_folder', '');           % beside it are in view
    if ~isempty(last) && isfolder(fileparts(last)), start = fileparts(last); end
    folder = uigetdir(start, 'Choose a test folder (test1_...) or a sample folder');
    if isequal(folder, 0), return; end
end
folder = char(folder);
if ~isfolder(folder), error('runDIC:folder', 'Not a folder: %s', folder); end
if usejava('desktop'), setpref('runDIC', 'last_folder', folder); end   % remembered across MATLAB sessions

% a test folder, or a sample folder holding test folders
[~, name] = fileparts(folder);
if ~isempty(regexp(name, '^test\d+_', 'once'))
    tests = {folder};
else
    d = dir(fullfile(folder, 'test*_*'));
    d = d([d.isdir] & ~cellfun(@isempty, regexp({d.name}, '^test\d+_', 'once')));
    tests = fullfile(folder, {d.name});
    if isempty(tests), error('runDIC:folder', 'No test folders (test1_..., test2_...) in %s', folder); end
end

if interactive && ~any(strcmp(varargin(1:2:end), 'steps'))
    s.steps = chooseSteps(tests, s.steps);
    if isempty(s.steps), util.log('Run cancelled at the step list, nothing run.'); return; end
end

try
    if any(strcmp(s.steps, 'info')) && usejava('desktop')
        askInfo(tests);                  % up front too
    end
    measuring = any(strcmp(s.steps, 'thickness'));
    if usejava('desktop') && (measuring || (any(strcmp(s.steps, 'export')) && s.ask_thickness))
        askThickness(tests, s, measuring);         % up front too, so the run doesn't stop later
    end
catch err
    if ~strcmp(err.identifier, 'runDIC:cancelled'), rethrow(err); end
    util.log('Run %s, nothing run.', err.message); return
end

s.n_tests = numel(tests);             % the checks offer to stop one test only when there are others
s.choices = containers.Map();          % the last tracking and export choices, offered again for the next test
tAll = tic;
done = false(size(tests));
for i = 1:numel(tests)
    t0 = tic;
    [~, tn] = fileparts(tests{i});
    util.log('\n=== %s (%d of %d)', tn, i, numel(tests));
    try
        runSteps(tests{i}, s);
        done(i) = true;
        util.log('=== %s done in %.1f min', tn, toc(t0) / 60);
    catch err
        if strcmp(err.identifier, 'runDIC:cancelled')      % the whole run, on purpose
            util.log('=== %s: run %s; %d test(s) not run', tn, err.message, numel(tests) - i);
            util.log('\nCancelled after %.1f min.', toc(tAll) / 60); return
        end
        if strcmp(err.identifier, 'runDIC:stopped')        % stopped on purpose at a check
            util.log('=== %s STOPPED: %s', tn, err.message); continue
        end
        util.log('=== %s FAILED: %s', tn, err.message);
        if isscalar(tests), rethrow(err); end
    end
end
util.log('\nFinished in %.1f min.', toc(tAll) / 60);
if usejava('desktop') && (any(strcmp(s.steps, 'report')) || (s.review_tracking && any(ismember({'track', 'artifacts'}, s.steps))))
    for i = find(done)
        viewTracking(tests{i}, s);   % the tracking report, one window per test
    end
end
end


function runSteps(testRoot, s)
need = @(f) assert(isfile(fullfile(testRoot, f)), 'runDIC:missing', 'Missing %s: run the earlier steps first', f);
needTips = @() assert(hasField(testRoot, 'tine_tips_px'), 'runDIC:missing', 'No tine tips yet: run the tines step first');
askTracking = any(strcmp(s.steps, 'track')) && s.ask_tracking && usejava('desktop');
skipTrack = false;                                   % Skip in the tracking window: the tracks saved before are kept
for st = s.steps
    if askTracking && ~ismember(st{1}, {'info', 'thickness', 'sort', 'tines', 'scale', 'widths'})   % after the tine tips, scale and spans
        [s.track_deformations, s.loading_only, skipTrack] = chooseTracking(testRoot, s);
        askTracking = false;
    end
    switch st{1}
        case 'thickness'
            if ~usejava('desktop'), util.log('-- thickness: needs the MATLAB desktop, skipped'); end
        case 'info'
            if ~usejava('desktop'), util.log('-- sample label and notes: needs the MATLAB desktop, skipped'); end
        case 'sort'
            util.log('-- sorting the images'); sortImages(testRoot, s);
        case 'tines'
            util.log('-- finding the tine tips'); findTineTips(testRoot, s.label, s);
        case 'scale'
            needTips();
            util.log('-- image scale'); calibrateScale(testRoot, s);
        case 'widths'
            need('scale_calibration.json');
            util.log('-- tine spans'); measureTineSpans(testRoot, s);
        case 'track'
            if skipTrack, util.log('-- tracking: skipped, the tracks saved before are kept'); continue; end
            need('scale_calibration.json');
            util.log('-- tracking'); trackProtocols(testRoot, s);
            trackingStats(testRoot, s);
        case 'artifacts'
            if skipTrack, util.log('-- tracking artifact check: skipped with the tracking'); continue; end
            util.log('-- tracking artifact check');
            if artifactCheck(testRoot, s)
                util.log('-- tracking again, with the new excluded regions'); trackProtocols(testRoot, s);
                util.log('-- checking again'); artifactCheck(testRoot, s, true);
                trackingStats(testRoot, s);
            end
        case 'export'
            if s.ask_export && usejava('desktop')                % once the tracking is done
                [s.export_deformations, s.export_formats, skip] = chooseExport(testRoot, s);
                if skip, util.log('-- export: skipped, the files saved before are kept'); continue; end
            end
            util.log('-- export'); exportStressStretch(testRoot, s);
            if any(strcmp(s.export_formats, 'xlsx')), exportAllData(testRoot, s); end
        case 'report'
            trackingStats(testRoot, s);     % the window opens when every test is done
        otherwise
            error('runDIC:step', 'Unknown step ''%s''', st{1});
    end
end
end


function steps = chooseSteps(tests, defaults)
% Ask which steps to run, noting which ones the (first) test has already done.
names = {'info', 'thickness', 'sort', 'tines', 'scale', 'widths', 'track', 'artifacts', 'export', 'report'};
what = {'sample label and notes', 'measure the thickness from the photos', 'sort the images into one folder per cycle', 'find the tine tips', 'image scale', ...
        'tine spans (specimen widths)', 'DIC tracking', 'check the tracking for local errors', ...
        'write the stress-stretch CSVs and the AllData file', 'open the tracking report (no tracking)'};
t = tests{1};
done = [~isempty(util.sampleInfo(t)), startsWith(thicknessSource(t), 'measure'), ~isempty(dir(fullfile(t, 'Images', 'Def11_*'))), ...
        hasField(t, 'tine_tips_px'), isfile(fullfile(t, 'scale_calibration.json')), ...
        hasField(t, 'widths_mm'), ~isempty(dir(fullfile(t, 'DIC_results', 'Def1*.mat'))), ...
        isfile(fullfile(t, 'DIC_results', 'artifact_check.png')), isfile(fullfile(t, 'stress_stretch', 'all_protocols.csv')), isfile(fullfile(t, 'DIC_results', 'tracking_report.csv'))];
items = cell(size(names));
for k = 1:numel(names)
    items{k} = sprintf('%-10s %s%s', names{k}, what{k}, repmat('   (done)', 1, done(k)));
end
if numel(tests) > 1, note = sprintf(' (done: as in %s)', util.shortName(t)); else, note = ''; end
[sel, ok] = listdlg('ListString', items, 'SelectionMode', 'multiple', 'InitialValue', find(ismember(names, defaults)), ...
    'Name', 'runDIC', 'PromptString', {sprintf('Steps to run on %d test(s)%s:', numel(tests), note), ...
    'Ctrl-click (Cmd-click on a Mac) to pick several; each step uses the earlier ones.'}, 'ListSize', util.ternary(isunix && ~ismac, [560 240], [460 190]), ...
    'OKString', 'Run');
if ~ok, steps = {}; else, steps = names(sel); end
end


function askThickness(tests, s, every)
% Ask for the wall thickness of the tests, typed (e.g. from the thickness workbook) or measured from the
% test's photos (measureThickness), and save the answers in the tests' test_config.json. With every (the
% thickness step runs), every test is listed, filled in with its saved thickness; otherwise only the tests
% without one. A blank box, or a saved thickness left as it is, changes nothing; a test without a
% thickness has NaN stress columns.
if every
    shown = tests;
else
    shown = tests(cellfun(@(t) isnan(util.thickness(t, s)), tests));
end
if isempty(shown), return; end
n = numel(shown);
Hsaved = nan(n, 1); srcSaved = repmat({''}, n, 1);
for i = 1:n
    [t, srcSaved{i}] = util.savedThickness(shown{i});
    if ~isempty(t), Hsaved(i) = t.H; end
end
W = 560; rowH = 66; Hfig = 136 + n * rowH;               % in pixels; buttons as in MATLAB's own dialogs
fig = figure('Name', 'runDIC: thickness for the stresses', 'NumberTitle', 'off', 'Color', 'w', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Position', [200 200 W Hfig], 'Resize', 'off', 'CloseRequestFcn', @(~, ~) cancelRun());
movegui(fig, 'center');
uicontrol(fig, 'Style', 'text', 'Position', [20 Hfig - 50 W - 40 36], 'FontSize', 10, ...
          'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', ...
          'Wall thickness of each test (mm): type it (e.g. from the thickness workbook), or measure it from the photos.');
edt = gobjects(n, 1); src = repmat({'saved'}, n, 1); sds = nan(n, 1); info = gobjects(n, 1); cancelled = false;
rev = gobjects(n, 1); notes = cell(n, 1);
for i = 1:n
    [sampleDir, name] = fileparts(shown{i});
    [~, sample] = fileparts(sampleDir);
    y = Hfig - 60 - (i - 1) * rowH - 24;                % the row's controls; its note on the lines below
    uicontrol(fig, 'Style', 'text', 'Position', [20 y + 3 170 16], 'FontSize', 10, ...
              'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', ...
              sprintf('%s, %s', strrep(sample, '_', ' '), regexp(name, '^test\d+', 'match', 'once')));
    edt(i) = uicontrol(fig, 'Style', 'edit', 'Position', [195 y 70 22], ...
                       'String', fmt(Hsaved(i)), 'Callback', @(~, ~) typed(i));
    [~, ph] = util.thicknessPhotos(shown{i});
    hm = util.button(fig, 'Measure from the photos', [275 y], @(~, ~) measure(i), ...
                     'Enable', util.ternary(isempty(ph), 'off', 'on'));
    p = get(hm, 'Position');
    rev(i) = util.button(fig, 'Revert', [p(1) + p(3) + 8, y], @(~, ~) revert(i), 'Enable', 'off');   % on once the box is yellow
    if isfinite(Hsaved(i))
        notes{i} = sprintf('saved: %s', util.ternary(isempty(srcSaved{i}), 'no source', srcSaved{i}));
    else
        notes{i} = util.ternary(isempty(ph), 'no photos found', sprintf('%d photo(s)', numel(ph)));
    end
    info(i) = uicontrol(fig, 'Style', 'text', 'Position', [195 y - 36 W - 215 32], ...
                        'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'ForegroundColor', [0.4 0.4 0.4], 'String', notes{i});
end
uicontrol(fig, 'Style', 'text', 'Position', [20 50 W - 40 30], ...
          'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'ForegroundColor', [0.4 0.4 0.4], 'String', ...
          'OK saves the yellow boxes; Revert puts a box back as it was. A blank box changes nothing.');
util.buttonRow(fig, {'OK', 'Cancel'}, {@(~, ~) uiresume(fig), @(~, ~) cancelRun()}, W - 20, 16);
uiwait(fig);
if cancelled || ~isvalid(fig)
    if isvalid(fig), delete(fig); end
    util.cancel('the thickness window');
end
answers = get(edt, 'String'); if ~iscell(answers), answers = {answers}; end
delete(fig);
for i = 1:n
    a = strtrim(answers{i});
    H = str2double(strrep(a, ',', '.'));
    [~, name] = fileparts(shown{i});
    test = regexp(name, '^test\d+', 'match', 'once');
    if isempty(a) || strcmp(a, fmt(Hsaved(i)))                          % as saved: its source is kept
        if isfinite(Hsaved(i))
            util.log('  %s: kept %.3f mm (%s)', test, Hsaved(i), srcSaved{i});
        else
            util.log('  %s: no thickness entered (stress columns NaN)', test);
        end
    elseif ~isfinite(H) || H <= 0 || H > 20
        util.log('  %s: "%s" is not a thickness in mm, skipped%s', test, a, ...
                 util.ternary(isfinite(Hsaved(i)), sprintf(' (kept %.3f mm)', Hsaved(i)), ' (stress columns NaN)'));
    else
        if strcmp(src{i}, 'saved'), src{i} = 'typed'; sds(i) = NaN; end    % edited without Enter
        util.setThickness(shown{i}, H, src{i}, sds(i));
        util.log('  %s: thickness %.3f mm (%s) saved in test_config.json%s', test, H, src{i}, ...
                 util.ternary(isfinite(Hsaved(i)), sprintf(' (was %.3f mm)', Hsaved(i)), ''));
    end
end

    function measure(i)
        r = measureThickness(shown{i}, s);
        if ~isvalid(fig), return; end
        figure(fig);                              % back to this window, measured or not
        if isempty(r), return; end
        set(edt(i), 'String', fmt(r.H), 'BackgroundColor', [1 0.96 0.75]);   % not saved yet
        set(rev(i), 'Enable', 'on');
        src{i} = r.source; sds(i) = r.sd;
        set(info(i), 'String', sprintf('measured %.3f +- %.3f mm (%s); OK saves it%s', r.H, r.sd, r.source, ...
            util.ternary(isfinite(Hsaved(i)), sprintf(' (saved before: %.3f mm)', Hsaved(i)), '')));
    end

    function typed(i)
        src{i} = 'typed'; sds(i) = NaN;
        set(edt(i), 'BackgroundColor', [1 0.96 0.75]);
        set(rev(i), 'Enable', 'on');
        set(info(i), 'String', 'typed; OK saves it');
    end

    function revert(i)
        src{i} = 'saved'; sds(i) = NaN;
        set(edt(i), 'String', fmt(Hsaved(i)), 'BackgroundColor', 'w');
        set(rev(i), 'Enable', 'off');
        set(info(i), 'String', notes{i});
    end

    function cancelRun()
        cancelled = true;
        uiresume(fig);
    end
end


function str = fmt(H)
% A thickness as shown in the thickness window ('' for none).
if isfinite(H), str = sprintf('%.3f', H); else, str = ''; end
end


function [defs, loadingOnly, skip] = chooseTracking(testRoot, s)
% Ask which deformations of the test to track and how much of each cycle, filled in with the answer for
% the previous test of the run, else the settings (track_deformations, loading_only). When the test has
% tracks, Skip leaves the tracking out (skip true) and keeps them.
c = s.choices;
if isKey(c, 'defs'), s.track_deformations = c('defs'); s.loading_only = c('loading_only'); end
W = 660; Hfig = 420;                                      % in pixels; buttons as in MATLAB's own dialogs
fig = figure('Name', 'runDIC: what to track', 'NumberTitle', 'off', 'Color', 'w', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Position', [200 200 W Hfig], 'Resize', 'off', 'CloseRequestFcn', @(~, ~) cancelRun());
movegui(fig, 'center');
uicontrol(fig, 'Style', 'text', 'Position', [20 Hfig - 36 W - 40 18], 'FontSize', 10, 'HorizontalAlignment', 'left', ...
          'BackgroundColor', 'w', 'String', sprintf('Deformations of %s to track:', util.shortName(testRoot)));
d = dir(fullfile(testRoot, 'DIC_results', 'Def*.mat'));
tracked = unique(cellfun(@(x) sscanf(x, 'Def%d'), {d.name}));
box = gobjects(19, 1);
heads = {'Def01-Def10: preconditioning', 'Def11-Def19: analyzed'};
for col = 1:2
    uicontrol(fig, 'Style', 'text', 'Position', [20 + (col - 1) * 320 Hfig - 60 300 16], 'FontWeight', 'bold', ...
              'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', heads{col});
end
for k = 1:19
    col = 1 + (k > 10); r = k - 10 * (col - 1);
    name = sprintf('Def%02d', k);
    try name = util.defFolder(testRoot, k); catch, end         % the sorted folder's name
    box(k) = uicontrol(fig, 'Style', 'checkbox', 'Position', [20 + (col - 1) * 320 Hfig - 62 - r * 20 300 18], ...
                       'BackgroundColor', 'w', 'Value', ismember(k, s.track_deformations), ...
                       'String', sprintf('%s%s', name, repmat('   (tracked)', 1, ismember(k, tracked))));
end
yb = Hfig - 62 - 11 * 20 - 10;
quick = {'All', 1:19; 'Def11-Def19', 11:19; 'None', []};
xq = 20;
for k = 1:3
    hq = util.button(fig, quick{k, 1}, [xq yb], @(~, ~) pick(quick{k, 2}));
    p = get(hq, 'Position'); xq = p(1) + p(3) + 8;
end
grp = uibuttongroup(fig, 'Units', 'pixels', 'Position', [20 yb - 74 W - 40 62], 'BackgroundColor', 'w', ...
                    'Title', 'Frames of each deformation');
rbLoad = uicontrol(grp, 'Style', 'radiobutton', 'Position', [10 24 W - 70 18], 'BackgroundColor', 'w', ...
                   'String', 'Loading branch only: from the start to the force peak (what the CSVs use)');
rbAll = uicontrol(grp, 'Style', 'radiobutton', 'Position', [10 4 W - 70 18], 'BackgroundColor', 'w', ...
                  'String', 'Loading and unloading: the whole cycle (more frames, so longer)');
grp.SelectedObject = util.ternary(s.loading_only, rbLoad, rbAll);
skipLabel = repmat({'Skip'}, 1, ~isempty(tracked));
util.buttonRow(fig, [{'Track'}, skipLabel, {'Cancel'}], [{@(~, ~) uiresume(fig)}, repmat({@(~, ~) skipStep()}, 1, numel(skipLabel)), ...
               {@(~, ~) cancelRun()}], W - 20, 14);
cancelled = false; skip = false;
uiwait(fig);
if cancelled || ~isvalid(fig)
    if isvalid(fig), delete(fig); end
    util.cancel('the tracking choice');
end
defs = find(cell2mat(get(box, 'Value')))';
loadingOnly = grp.SelectedObject == rbLoad;
delete(fig);
if ~skip
    if isempty(defs), util.cancel('the tracking choice (no deformation picked)'); end
    c('defs') = defs; c('loading_only') = loadingOnly; %#ok<NASGU> a handle: kept for the next test
    util.log('-- tracking %s, %s', util.defList(defs), util.ternary(loadingOnly, 'loading branches only', 'whole cycles'));
end

    function pick(k)
        for j = 1:19, box(j).Value = ismember(j, k); end
    end

    function cancelRun()
        cancelled = true;
        uiresume(fig);
    end

    function skipStep()
        skip = true;
        uiresume(fig);
    end
end


function [defs, formats, skip] = chooseExport(testRoot, s)
% Ask which deformations of the test to export and which files to write, filled in with the answer for
% the previous test of the run, else the settings (export_deformations, export_formats). Deformations not
% tracked are grayed out. When files were exported before, Skip leaves the export out (skip true).
c = s.choices;
if isKey(c, 'export_defs'), s.export_deformations = c('export_defs'); s.export_formats = c('export_formats'); end
d = dir(fullfile(testRoot, 'DIC_results', 'Def*.mat'));
avail = unique(cellfun(@(x) sscanf(x, 'Def%d'), {d.name}));
fmt = {'csv', 'One CSV per protocol (DefNN_<name>.csv)'; 'all_csv', 'All protocols in one CSV (all_protocols.csv)'; ...
       'png', 'Stress-stretch plot, PNG'; 'pdf', 'Stress-stretch plot, PDF'; 'fig', 'Stress-stretch plot, MATLAB figure (.fig)'; ...
       'xlsx', 'AllData file, CellScale''s layout (.xlsx)'};
W = 660; Hfig = 500;                                      % in pixels; buttons as in MATLAB's own dialogs
fig = figure('Name', 'runDIC: what to export', 'NumberTitle', 'off', 'Color', 'w', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Position', [200 200 W Hfig], 'Resize', 'off', 'CloseRequestFcn', @(~, ~) cancelRun());
movegui(fig, 'center');
uicontrol(fig, 'Style', 'text', 'Position', [20 Hfig - 36 W - 40 18], 'FontSize', 10, 'HorizontalAlignment', 'left', ...
          'BackgroundColor', 'w', 'String', sprintf('Deformations of %s to export (grayed out: not tracked):', util.shortName(testRoot)));
box = gobjects(19, 1);
heads = {'Def01-Def10: preconditioning', 'Def11-Def19: analyzed'};
for col = 1:2
    uicontrol(fig, 'Style', 'text', 'Position', [20 + (col - 1) * 320 Hfig - 60 300 16], 'FontWeight', 'bold', ...
              'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', heads{col});
end
for k = 1:19
    col = 1 + (k > 10); r = k - 10 * (col - 1);
    name = sprintf('Def%02d', k);
    try name = util.defFolder(testRoot, k); catch, end         % the sorted folder's name
    box(k) = uicontrol(fig, 'Style', 'checkbox', 'Position', [20 + (col - 1) * 320 Hfig - 62 - r * 20 300 18], ...
                       'BackgroundColor', 'w', 'String', name, 'Value', ismember(k, s.export_deformations) && ismember(k, avail), ...
                       'Enable', util.ternary(ismember(k, avail), 'on', 'off'));
end
yb = Hfig - 62 - 11 * 20 - 10;
quick = {'All', 1:19; 'Def11-Def19', 11:19; 'None', []};
xq = 20;
for k = 1:3
    hq = util.button(fig, quick{k, 1}, [xq yb], @(~, ~) pick(quick{k, 2}));
    p = get(hq, 'Position'); xq = p(1) + p(3) + 8;
end
nf = size(fmt, 1);
grp = uipanel(fig, 'Units', 'pixels', 'Position', [20 yb - 30 - nf * 20 W - 40 nf * 20 + 22], 'BackgroundColor', 'w', ...
              'Title', 'Files (README.txt is always written)');
fbox = gobjects(nf, 1);
for k = 1:nf
    fbox(k) = uicontrol(grp, 'Style', 'checkbox', 'Position', [10 (nf - k) * 20 + 4 W - 70 18], 'BackgroundColor', 'w', ...
                        'String', fmt{k, 2}, 'Value', ismember(fmt{k, 1}, s.export_formats));
end
old = [dir(fullfile(testRoot, 'stress_stretch', '*.csv')); dir(fullfile(testRoot, 'stress_stretch', '*.xlsx'))];
skipLabel = repmat({'Skip'}, 1, ~isempty(old));
util.buttonRow(fig, [{'Export'}, skipLabel, {'Cancel'}], [{@(~, ~) uiresume(fig)}, repmat({@(~, ~) skipStep()}, 1, numel(skipLabel)), ...
               {@(~, ~) cancelRun()}], W - 20, 14);
cancelled = false; skip = false;
uiwait(fig);
if cancelled || ~isvalid(fig)
    if isvalid(fig), delete(fig); end
    util.cancel('the export choice');
end
defs = find(cell2mat(get(box, 'Value')))';
formats = fmt(logical(cell2mat(get(fbox, 'Value'))), 1)';
delete(fig);
if ~skip
    c('export_defs') = defs; c('export_formats') = formats; %#ok<NASGU> a handle: kept for the next test
    if isempty(defs) || isempty(formats)
        util.log('-- export: nothing picked, only README.txt');
    else
        util.log('-- export %s: %s', util.defList(defs), strjoin(formats, ', '));
    end
end

    function pick(k)
        for j = 1:19, box(j).Value = ismember(j, k) && ismember(j, avail); end
    end

    function cancelRun()
        cancelled = true;
        uiresume(fig);
    end

    function skipStep()
        skip = true;
        uiresume(fig);
    end
end


function src = thicknessSource(t)
% How the test's saved thickness was obtained ('' if none).
[~, src] = util.savedThickness(t);
end


function askInfo(tests)
% Ask for the label and notes of every test (filled in with the saved ones) and save them in each
% test's test_config.json.
prompts = cell(1, 2 * numel(tests)); defaults = prompts;
for i = 1:numel(tests)
    [sampleDir, name] = fileparts(tests{i});
    [~, sample] = fileparts(sampleDir);
    who = sprintf('%s, %s', strrep(sample, '_', ' '), regexp(name, '^test\d+', 'match', 'once'));
    [label, notes] = util.sampleInfo(tests{i});
    prompts{2 * i - 1} = sprintf('%s: label (for example healthy, or slightly calcified)', who);
    prompts{2 * i} = sprintf('%s: notes, one per line (tears, re-mounting, anything unusual)', who);
    defaults{2 * i - 1} = label; defaults{2 * i} = char(notes);
end
answers = inputdlg(prompts, 'runDIC: sample label and notes', repmat([1 80; 4 80], numel(tests), 1), defaults, ...
                   struct('Resize', 'on'));
if isempty(answers), util.cancel('the label and notes window'); end
for i = 1:numel(tests)
    label = strtrim(answers{2 * i - 1});
    notes = strtrim(cellstr(answers{2 * i}))';
    notes = notes(~cellfun(@isempty, notes));
    util.setSampleInfo(tests{i}, label, notes);
    util.log('  %s: label "%s", %d note(s) -> test_config.json', util.shortName(tests{i}), label, numel(notes));
end
end


function ok = hasField(t, field)
% True when the test's test_config.json has `field`.
ok = false;
f = fullfile(t, 'test_config.json');
if isfile(f), c = util.readJSON(f); ok = isfield(c, field); end
end
