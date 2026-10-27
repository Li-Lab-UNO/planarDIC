function runDIC(folder, varargin)
% runDIC  Stress-stretch data from a CellScale BioTester planar biaxial test.
%
%   runDIC              asks for a folder (a test folder test1_..., or a sample folder holding tests),
%                       then which steps to run
%   runDIC(folder)      runs that test, or every test in that sample folder
%   runDIC(folder, 'steps', {'track', 'export'})   only some steps
%   runDIC(folder, 'thickness_mm', 1.25)           thickness for the stresses (mm)
%
% Settings are in DIC_settings.m. A test folder holds the CellScale output: *Data.csv, Coords1.csv,
% Logs/*_PhaseTiming.csv and the images in Images/. The steps:
%   info    asks for the sample label (e.g. healthy) and notes -> test_config.json
%   thickness  measures the wall thickness from the test's photos -> test_config.json; if
%           the test already has one, asks whether to keep it or measure again
%   sort    images into one folder per load-unload cycle (Images/Def01_... to Def19_...), then checks that
%           every cycle has all its images, Data.csv rows and Coords1.csv rows -> Images/sort_summary.csv
%   tines   finds the 20 tine tips in the first image of Def11 -> test_config.json
%           check DIC_results/config_tine_tips.png: every red circle should sit on a tine tip
%   scale   image scale from the tine motion in Def01 and Def19 -> scale_calibration.json
%           check DIC_results/scale_calibration.png
%   widths  tine spans (specimen widths for the stresses) -> test_config.json
%   branches  shows the loading branches of Def11-Def19 on the rake forces; wrong ones can be clicked
%           -> loading_branches.json, DIC_results/loading_branches.png
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
    if isempty(s.steps), return; end
end

if any(strcmp(s.steps, 'info')) && usejava('desktop')
    askInfo(tests);                  % up front too
end
if any(strcmp(s.steps, 'export')) && s.ask_thickness && usejava('desktop')
    askThickness(tests, s, any(strcmp(s.steps, 'thickness')));          % up front, so the run doesn't stop later
end

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
for st = s.steps
    switch st{1}
        case 'thickness'
            util.log('-- thickness from the photos'); thicknessStep(testRoot, s);
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
        case 'branches'
            util.log('-- loading branches'); reviewBranches(testRoot, s);
        case 'track'
            need('scale_calibration.json');
            util.log('-- tracking (the long step)'); trackProtocols(testRoot, s);
            trackingStats(testRoot, s);
        case 'artifacts'
            util.log('-- tracking artifact check');
            if artifactCheck(testRoot, s)
                util.log('-- tracking again, with the new excluded regions'); trackProtocols(testRoot, s);
                util.log('-- checking again'); artifactCheck(testRoot, s, true);
                trackingStats(testRoot, s);
            end
        case 'export'
            util.log('-- stress-stretch CSVs'); exportStressStretch(testRoot, s); exportAllData(testRoot, s);
        case 'report'
            trackingStats(testRoot, s);     % the window opens when every test is done
        otherwise
            error('runDIC:step', 'Unknown step ''%s''', st{1});
    end
end
end


function steps = chooseSteps(tests, defaults)
% Ask which steps to run, noting which ones the (first) test has already done.
names = {'info', 'thickness', 'sort', 'tines', 'scale', 'widths', 'branches', 'track', 'artifacts', 'export', 'report'};
what = {'sample label and notes', 'measure the thickness from the photos', 'sort the images into one folder per cycle', 'find the tine tips', 'image scale', ...
        'tine spans (specimen widths)', 'check the loading branches', 'DIC tracking (the long step)', 'check the tracking for local errors', ...
        'write the stress-stretch CSVs and the AllData file', 'open the tracking report (no tracking)'};
t = tests{1};
done = [~isempty(util.sampleInfo(t)), startsWith(thicknessSource(t), 'measure'), ~isempty(dir(fullfile(t, 'Images', 'Def11_*'))), ...
        hasField(t, 'tine_tips_px'), isfile(fullfile(t, 'scale_calibration.json')), ...
        hasField(t, 'widths_mm'), isfile(fullfile(t, 'loading_branches.json')), ~isempty(dir(fullfile(t, 'DIC_results', 'Def1*.mat'))), ...
        isfile(fullfile(t, 'DIC_results', 'artifact_check.png')), isfile(fullfile(t, 'stress_stretch', 'all_protocols.csv')), isfile(fullfile(t, 'DIC_results', 'tracking_report.csv'))];
items = cell(size(names));
for k = 1:numel(names)
    items{k} = sprintf('%-10s %s%s', names{k}, what{k}, repmat('   (done)', 1, done(k)));
end
if numel(tests) > 1, note = sprintf(' (done: as in %s)', util.shortName(t)); else, note = ''; end
[sel, ok] = listdlg('ListString', items, 'SelectionMode', 'multiple', 'InitialValue', find(ismember(names, defaults)), ...
    'Name', 'runDIC', 'PromptString', {sprintf('Steps to run on %d test(s)%s:', numel(tests), note), ...
    'Ctrl-click (Cmd-click on a Mac) to pick several; each step uses the earlier ones.'}, 'ListSize', [460 190], ...
    'OKString', 'Run');
if ~ok, steps = {}; else, steps = names(sel); end
end


function askThickness(tests, s, measuring)
% Ask for the thickness of every test that has none, typed or measured from the test's photos
% (measureThickness), and save the answers in the tests' test_config.json. Blank answers are
% skipped (those stress columns stay NaN). With measuring (the thickness step runs), tests with photos
% are left to that step.
missing = {};
for i = 1:numel(tests)
    [~, ph] = util.thicknessPhotos(tests{i});
    if isnan(util.thickness(tests{i}, s)) && ~(measuring && ~isempty(ph)), missing{end + 1} = tests{i}; end %#ok<AGROW>
end
if isempty(missing), return; end
n = numel(missing); rowH = 0.8 / max(n, 3);
fig = figure('Name', 'runDIC: thickness for the stresses', 'NumberTitle', 'off', 'Color', 'w', 'MenuBar', 'none', ...
             'ToolBar', 'none', 'Position', [200 200 820 120 + 60 * n], 'CloseRequestFcn', @(f, ~) uiresume(f));
movegui(fig, 'center');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.03 0.86 0.94 0.1], 'FontSize', 12, ...
          'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', ...
          'Wall thickness of each test (mm): type it, or measure it from the photos. Leave blank to skip.');
edt = gobjects(n, 1); src = repmat({'typed'}, n, 1); sds = nan(n, 1); info = gobjects(n, 1);
for i = 1:n
    [sampleDir, name] = fileparts(missing{i});
    [~, sample] = fileparts(sampleDir);
    y = 0.84 - i * rowH;
    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.03 y 0.3 rowH * 0.6], 'FontSize', 12, ...
              'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'String', ...
              sprintf('%s, %s', strrep(sample, '_', ' '), regexp(name, '^test\d+', 'match', 'once')));
    edt(i) = uicontrol(fig, 'Style', 'edit', 'Units', 'normalized', 'Position', [0.33 y + rowH * 0.1 0.12 rowH * 0.55], ...
                       'FontSize', 12, 'Callback', @(~, ~) typed(i));
    [~, ph] = util.thicknessPhotos(missing{i});
    uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.47 y + rowH * 0.1 0.24 rowH * 0.55], ...
              'FontSize', 12, 'String', 'Measure from the photos', 'Enable', util.ternary(isempty(ph), 'off', 'on'), ...
              'Callback', @(~, ~) measure(i));
    info(i) = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.73 y 0.25 rowH * 0.6], 'FontSize', 10, ...
                        'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'ForegroundColor', [0.4 0.4 0.4], ...
                        'String', util.ternary(isempty(ph), 'no photos found', sprintf('%d photo(s)', numel(ph))));
end
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.62 0.03 0.16 0.1], 'String', 'OK', ...
          'FontSize', 12, 'Callback', @(~, ~) uiresume(fig));
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.8 0.03 0.16 0.1], 'String', 'Skip all', ...
          'FontSize', 12, 'Callback', @(~, ~) skipAll());
uiwait(fig);
if ~isvalid(fig), return; end
answers = get(edt, 'String'); if ~iscell(answers), answers = {answers}; end
delete(fig);
for i = 1:n
    H = str2double(strrep(strtrim(answers{i}), ',', '.'));
    [~, name] = fileparts(missing{i});
    test = regexp(name, '^test\d+', 'match', 'once');
    if isempty(strtrim(answers{i}))
        util.log('  %s: no thickness entered (stress columns NaN)', test);
    elseif ~isfinite(H) || H <= 0 || H > 20
        util.log('  %s: "%s" is not a thickness in mm, skipped (stress columns NaN)', test, answers{i});
    else
        util.setThickness(missing{i}, H, src{i}, sds(i));
        util.log('  %s: thickness %.3f mm (%s) saved in test_config.json', test, H, src{i});
    end
end

    function measure(i)
        r = measureThickness(missing{i}, s);
        if isempty(r) || ~isvalid(fig), return; end
        set(edt(i), 'String', sprintf('%.3f', r.H));
        src{i} = r.source; sds(i) = r.sd;
        set(info(i), 'String', sprintf('%.3f +- %.3f mm, %s', r.H, r.sd, r.source));
        figure(fig);
    end

    function typed(i)
        src{i} = 'typed'; sds(i) = NaN;
        set(info(i), 'String', 'typed');
    end

    function skipAll()
        set(edt, 'String', '');
        uiresume(fig);
    end
end


function thicknessStep(testRoot, s)
% Measure the thickness from the test's photos and save it in the test's test_config.json. When the
% test already has one, ask first whether to keep it (the default) or measure again.
[~, name] = fileparts(testRoot);
test = regexp(name, '^test\d+', 'match', 'once');
[~, photos] = util.thicknessPhotos(testRoot);
if isempty(photos), util.log('  no thickness photos (Thickness*\\test %s), skipped', test(5:end)); return; end
if ~usejava('desktop'), util.log('  measuring needs the MATLAB desktop, skipped'); return; end
s.thickness_mm = [];                          % the saved thickness, not one given to runDIC
[Hold, srcOld] = util.thickness(testRoot, s);
if isfinite(Hold)
    q = questdlg({sprintf('%s already has a thickness: %.3f mm', test, Hold), sprintf('(%s).', srcOld), '', ...
                  'Keep it, or measure it again from the photos?'}, 'runDIC: thickness', ...
                 'Keep it and proceed', 'Measure again', 'Keep it and proceed');
    if ~strcmp(q, 'Measure again'), util.log('  kept %.3f mm (%s)', Hold, srcOld); return; end
end
r = measureThickness(testRoot, s);
if isempty(r)
    if isfinite(Hold), util.log('  nothing measured, kept %.3f mm (%s)', Hold, srcOld); return; end
    a = inputdlg(sprintf('No thickness measured for %s. Type it (mm), or leave blank to skip:', test), ...
                 'runDIC: thickness', [1 60]);
    H = str2double(strrep(strtrim([a{:} '']), ',', '.'));
    if isfinite(H) && H > 0 && H <= 20
        util.setThickness(testRoot, H, 'typed');
        util.log('  %s: thickness %.3f mm (typed) saved in test_config.json', test, H);
    else
        util.log('  %s: no thickness (stress columns NaN)', test);
    end
    return
end
util.setThickness(testRoot, r.H, r.source, r.sd);
util.log('  %s: thickness %.3f +- %.3f mm (%s) saved in test_config.json%s', test, r.H, r.sd, r.source, ...
         util.ternary(isfinite(Hold), sprintf(' (was %.3f mm)', Hold), ''));
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
if isempty(answers), util.log('Sample label and notes: cancelled, nothing changed'); return; end
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
