function viewTracking(testRoot, settings)
% viewTracking  Tracking report of a test: the analyzed deformations (Def11-Def19) side by side, played
%               back together.
%
%   viewTracking            asks for a test folder (opens where the last runDIC folder is)
%   viewTracking(folder)    that test folder
%
% runDIC opens it at the end of a run that tracked (review_tracking in DIC_settings.m). Each tile shows
% one deformation: the tracked points over the images (green: tracked, red: failed in that frame and
% filled from the neighbors), every 4th grid row and column, and the stretches, valid points and match
% quality (ZNCC) of the frame shown. One slider moves all of them by percent of their loading branch.
% Tiles marked Problem (red) or Check (amber) fall below the limits in DIC_settings.m; the window opens
% at the worst frame of the worst one. "Open large" shows the selected deformation in a big viewer.
% Images are read as they are shown, so it needs little memory; playback keeps time on slow computers
% by showing fewer frames.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
if nargin < 2 || isempty(settings), settings = DIC_settings(); end
if nargin < 1 || isempty(testRoot)
    start = pwd;
    last = getpref('runDIC', 'last_folder', '');
    if ~isempty(last) && isfolder(last), start = last; end
    testRoot = uigetdir(start, 'Choose a test folder (test1_...)');
    if isequal(testRoot, 0), return; end
end
testRoot = char(testRoot);
R = trackingStats(testRoot, settings, true);
R = R([R.analyzed]);                             % the window shows the analyzed ones (Def11-Def19)
cfg = util.readJSON(fullfile(testRoot, 'test_config.json'));
n = numel(R);
nExcl = 0;
if isfield(cfg, 'dic') && isfield(cfg.dic, 'exclude_px'), nExcl = numel(cfg.dic.exclude_px) / 4; end

LINE_EVERY = 4;
TILE_PX = 360;                                   % displayed tile size (image px): images are subsampled to it
speeds = [1 2 5 10];
COL = struct('good', [0.15 0.55 0.15], 'check', [0.85 0.55 0], 'problem', [0.85 0.1 0.1]);

% ---- start: the worst frame of the worst deformation, or the peak when all are good
sev = cellfun(@(s) find(strcmp(s, {'good', 'check', 'problem'})), {R.status});
[~, sel] = sortrows([-sev(:), [R.valid_worst]'], [1 2]); sel = sel(1);
if sev(sel) == 1, f = 1; else, f = R(sel).worst_fraction; end
playing = false;

% ---- window
fig = figure('Name', ['runDIC: tracking report, ' util.shortName(testRoot)], 'NumberTitle', 'off', 'Color', 'w', ...
             'WindowState', 'maximized', 'MenuBar', 'none', 'ToolBar', 'none', 'CloseRequestFcn', @onClose);
[msg, sub, mcol] = verdictText();
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.955 0.8 0.035], 'String', msg, ...
          'FontSize', 13, 'FontWeight', 'bold', 'ForegroundColor', mcol, 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.01 0.925 0.8 0.03], 'String', sub, ...
          'FontSize', 11, 'ForegroundColor', [0.35 0.35 0.35], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
chkGrid = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.84 0.95 0.07 0.03], 'String', 'Grid', ...
                    'Value', 1, 'FontSize', 11, 'BackgroundColor', 'w', 'Callback', @(~, ~) drawAll(true));
chkPts = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.91 0.95 0.08 0.03], 'String', 'Points', ...
                   'Value', 1, 'FontSize', 11, 'BackgroundColor', 'w', 'Callback', @(~, ~) drawAll(true));

tiles = cell(n, 1);
nc = ceil(sqrt(n)); nr = ceil(n / nc);
for i = 1:n
    tiles{i} = makeTile(i, mod(i - 1, nc), floor((i - 1) / nc));
end

btnPlay = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.01 0.025 0.06 0.04], ...
                    'String', 'Play', 'FontSize', 12, 'Callback', @(~, ~) onPlay());
popSpeed = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.075 0.025 0.07 0.04], ...
                     'String', compose('%dx', speeds), 'Value', 2, 'FontSize', 11);
slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.16 0.035 0.6 0.03], ...
                   'Min', 0, 'Max', 1, 'Value', f, 'SliderStep', [0.01 0.03], 'Callback', @(~, ~) onSlider());
addlistener(slider, 'ContinuousValueChange', @(~, ~) onSlider());
axMark = axes(fig, 'Units', 'normalized', 'Position', [0.172 0.012 0.576 0.02], 'XLim', [0 1], 'YLim', [0 1], ...
              'Visible', 'off'); hold(axMark, 'on');
for i = find(sev > 1)
    plot(axMark, R(i).worst_fraction, 0.5, '^', 'MarkerSize', 8, 'MarkerFaceColor', COL.(R(i).status), ...
         'MarkerEdgeColor', 'none');
end
txtPct = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.77 0.028 0.12 0.03], 'FontSize', 11, ...
                   'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.9 0.025 0.09 0.04], ...
          'String', 'Open large', 'FontSize', 12, 'Callback', @(~, ~) onLarge());
drawAll(true);


    function T = makeTile(i, col, row)
        % One deformation: axes, image (the tine region, subsampled), grid lines, points, readout.
        w = 0.98 / nc; h = 0.79 / nr;
        ax = axes(fig, 'Units', 'normalized', 'Position', [0.01 + col * w + 0.01, 0.095 + (nr - 1 - row) * h + 0.05, w - 0.02, h - 0.085]);
        T = struct('ax', ax);
        T.trk = load(fullfile(testRoot, 'DIC_results', [R(i).name '.mat']), 'posx', 'posy', 'valid', 'zncc', ...
                     'frames', 'times', 'X0', 'shape');
        T.imgDir = fullfile(testRoot, 'Images', R(i).name);
        info = imfinfo(fullfile(T.imgDir, [T.trk.frames{1} '.jpg']));
        z = cfg.tine_tips_px.z + 1; t = cfg.tine_tips_px.theta + 1;
        m = 0.35 * max(diff(z), diff(t));
        T.cols = [max(1, round(z(1) - m)), min(info.Width, round(z(2) + m))];
        T.rows = [max(1, round(t(1) - m)), min(info.Height, round(t(2) + m))];
        T.ds = max(1, floor(max(diff(T.cols), diff(T.rows)) / TILE_PX));
        T.hImg = image(ax, zeros(2, 2, 'uint8'), 'CDataMapping', 'scaled', 'ButtonDownFcn', @(~, ~) onSelect(i));
        colormap(ax, gray); set(ax, 'CLim', [0 255]); axis(ax, 'image'); hold(ax, 'on');
        toTile = @(v, c0) (v + 1 - c0) / T.ds + 1;            % native 0-based px -> tile px
        T.toX = @(x) toTile(x, T.cols(1)); T.toY = @(y) toTile(y, T.rows(1));
        plot(ax, T.toX(cfg.tine_tips_px.z([1 2 2 1 1])), T.toY(cfg.tine_tips_px.theta([1 1 2 2 1])), '--', ...
             'Color', [1 0.75 0.1], 'PickableParts', 'none');
        T.hGrid = plot(ax, NaN, NaN, '-', 'Color', [0.3 0.9 1 0.6], 'LineWidth', 0.6, 'PickableParts', 'none');
        T.hGood = plot(ax, NaN, NaN, '.', 'Color', [0.2 0.9 0.2], 'MarkerSize', 4, 'PickableParts', 'none');
        T.hBad = plot(ax, NaN, NaN, 'x', 'Color', [1 0.2 0.2], 'MarkerSize', 6, 'LineWidth', 1.2, 'PickableParts', 'none');
        set(ax, 'XTick', [], 'YTick', [], 'Box', 'on');
        T.lastK = 0;
        [cx, cy] = ndgrid(1:T.trk.shape(2), 1:T.trk.shape(1));            % column, row of each point
        T.sub = mod(cx(:), 2) == 1 & mod(cy(:), 2) == 1;
    end

    function drawAll(force)
        if ~isvalid(fig), return; end
        for q = 1:n, drawTile(q, force); end
        set(txtPct, 'String', sprintf('%.0f%% of loading', 100 * f));
        drawnow
    end

    function drawTile(i, force)
        T = tiles{i}; trk = T.trk; nT = numel(trk.times);
        k = round(f * (nT - 1)) + 1;
        if k ~= T.lastK
            set(T.hImg, 'CData', readCrop(fullfile(T.imgDir, [trk.frames{k} '.jpg']), T.rows, T.cols, T.ds));
            if T.lastK == 0, axis(T.ax, 'image'); end
            tiles{i}.lastK = k;
        elseif ~force
            return
        end
        x = trk.posx(k, :)'; y = trk.posy(k, :)'; v = trk.valid(k, :)';
        if get(chkPts, 'Value')                       % every 2nd row and column of the tracked points (the
            g = v & T.sub;                             % tiles are small); every failed point
            set(T.hGood, 'XData', T.toX(x(g)), 'YData', T.toY(y(g)));
            set(T.hBad, 'XData', T.toX(x(~v)), 'YData', T.toY(y(~v)));
        else
            set([T.hGood T.hBad], 'XData', NaN, 'YData', NaN);
        end
        if get(chkGrid, 'Value')
            [gx, gy] = trackGridLines(x, y, trk.shape, LINE_EVERY);
            set(T.hGrid, 'XData', T.toX(gx), 'YData', T.toY(gy));
        else
            set(T.hGrid, 'XData', NaN, 'YData', NaN);
        end
        F = dicCore.affineFit(trk.X0(v, :), [x(v), y(v)]); C = F' * F;
        xlabel(T.ax, {sprintf('frame %d/%d   valid %.1f%%   ZNCC %.3f', k, nT, 100 * mean(v), ...
                              median(trk.zncc(k, :), 'omitnan')), ...
                      sprintf('lambda_z %.3f   lambda_theta %.3f', sqrt(C(1, 1)), sqrt(C(2, 2)))}, ...
               'Interpreter', 'none', 'FontSize', 9, 'FontName', 'FixedWidth', 'Color', [0.2 0.2 0.2]);
        c = COL.(R(i).status);
        lab = sprintf('%s   %s', strrep(R(i).name, '_', ' '), upper(R(i).status));
        if i == sel, lab = ['> ' lab]; end
        title(T.ax, lab, 'Color', c, 'FontSize', 12, 'FontWeight', util.ternary(i == sel, 'bold', 'normal'), 'Interpreter', 'none');
        set(T.ax, 'XColor', c, 'YColor', c, 'LineWidth', util.ternary(strcmp(R(i).status, 'good'), 0.5, 3));
        T.ax.XLabel.Color = [0.2 0.2 0.2];
    end

    function onSlider()
        f = get(slider, 'Value');
        drawAll(false);
    end

    function onSelect(i)
        old = sel; sel = i;
        drawTile(old, true); drawTile(sel, true); drawnow
    end

    function onPlay()
        playing = ~playing;
        set(btnPlay, 'String', util.ternary(playing, 'Pause', 'Play'));
        if ~playing, return; end
        if f >= 1, f = 0; end
        dur = max(arrayfun(@(c) c{1}.trk.times(end) - c{1}.trk.times(1), tiles));   % longest branch (s)
        f0 = f; t0 = tic;
        while playing && isvalid(fig) && f < 1
            f = min(1, f0 + toc(t0) * speeds(get(popSpeed, 'Value')) / dur);    % paced by the clock
            set(slider, 'Value', f);
            drawAll(false);
        end
        if isvalid(fig), playing = false; set(btnPlay, 'String', 'Play'); end
    end

    function onLarge()
        T = tiles{sel};
        viewDeformation(testRoot, R(sel).name, round(f * (numel(T.trk.times) - 1)) + 1);
    end

    function onClose(~, ~)
        playing = false;
        delete(fig);
    end

    function [msg, sub, c] = verdictText()
        nb = sum(sev == 3); nck = sum(sev == 2);
        if nExcl == 0, art = 'Artifact check: nothing excluded.';
        else, art = sprintf('Artifact check: %d region(s) excluded and tracked again.', nExcl); end
        if nb + nck == 0
            msg = sprintf('All %d deformations tracked well. %s', n, art); c = COL.good;
            sub = 'Press Play to watch them; the slider moves all of them through their loading branch.';
        else
            p = {};
            if nb, p{end + 1} = sprintf('%d %s a tracking problem', nb, util.ternary(nb == 1, 'deformation has', 'deformations have')); end
            if nck && nb, p{end + 1} = sprintf('%d %s a look', nck, util.ternary(nck == 1, 'needs', 'need'));
            elseif nck, p{end + 1} = sprintf('%d %s a look', nck, util.ternary(nck == 1, 'deformation needs', 'deformations need')); end
            msg = sprintf('%s. %s', strjoin(p, ', '), art);
            c = util.ternary(nb > 0, COL.problem, COL.check);
            sub = sprintf(['Opened at the worst frame of %s. Play it, and show this window to your supervisor ' ...
                           'before using its results.'], strtok(R(sel).name, '_'));
        end
    end
end


function im = readCrop(path, rows, cols, ds)
% Gray image (uint8) of the rows/cols range of a jpg, every ds-th pixel.
im = imread(path);
im = im(rows(1):ds:rows(2), cols(1):ds:cols(2), :);
if size(im, 3) == 3
    im = uint8(0.299 * double(im(:, :, 1)) + 0.587 * double(im(:, :, 2)) + 0.114 * double(im(:, :, 3)));
end
end
