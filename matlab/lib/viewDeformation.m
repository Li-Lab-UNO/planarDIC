function viewDeformation(testRoot, startName, startFrame)
% Play back the DIC tracking of one deformation at a time, large: the tracked points over the images,
% frame by frame (opened by "Open large" in the tracking report, viewTracking). startName and
% startFrame choose the deformation and frame it opens at.
%
% Green dots are points tracked in that frame; red crosses are points that failed there and were filled
% from their neighbors. The grid lines (every 4th row and column) show how the tissue deforms. The
% stretches are those of the export (affine fit of all valid points, relative to the first image of
% the reference deformation). It opens zoomed on the tine region; zoom and pan with the axes toolbar.
% Drag the slider, or press Play.

testRoot = char(testRoot);
d = dir(fullfile(testRoot, 'DIC_results', 'Def*.mat'));
d = d(cellfun(@isempty, regexp({d.name}, '_fields\.mat$', 'once')));
if isempty(d), error('runDIC:missing', 'No tracks in %s: run the track step first', fullfile(testRoot, 'DIC_results')); end
names = regexprep({d.name}, '\.mat$', '');
cfg = struct();
if isfile(fullfile(testRoot, 'test_config.json')), cfg = util.readJSON(fullfile(testRoot, 'test_config.json')); end

LINE_EVERY = 4;                                     % grid lines drawn every 4th row and column
speeds = [1 2 5 10];                                % playback, times real time
trk = []; k = 1; playing = false; imgDir = '';

% ---- window
fig = figure('Name', ['runDIC: tracking, ' util.shortName(testRoot)], 'NumberTitle', 'off', 'Color', 'w', ...
             'WindowState', 'maximized', 'CloseRequestFcn', @onClose);
ax = axes(fig, 'Units', 'normalized', 'Position', [0.03 0.1 0.72 0.86]);
hImg = imagesc(ax, zeros(10)); colormap(ax, gray); axis(ax, 'image', 'off'); hold(ax, 'on');
hTips = plot(ax, NaN, NaN, '--', 'Color', [1 0.75 0.1], 'LineWidth', 1);
hGrid = plot(ax, NaN, NaN, '-', 'Color', [0.3 0.9 1 0.6], 'LineWidth', 0.6);
hGood = plot(ax, NaN, NaN, '.', 'Color', [0.2 0.9 0.2], 'MarkerSize', 5);
hBad = plot(ax, NaN, NaN, 'x', 'Color', [1 0.2 0.2], 'MarkerSize', 7, 'LineWidth', 1.5);
slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.03 0.04 0.72 0.03], ...
                   'Min', 1, 'Max', 2, 'Value', 1, 'Callback', @(~, ~) onSlider());
addlistener(slider, 'ContinuousValueChange', @(~, ~) onSlider());

x0 = 0.78; w = 0.2;
label(0.93, 'Deformation');
popDef = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [x0 0.89 w 0.04], ...
                   'String', names, 'FontSize', 11, 'Callback', @(~, ~) loadDef());
btnPlay = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [x0 0.82 w * 0.48 0.05], ...
                    'String', 'Play', 'FontSize', 12, 'Callback', @(~, ~) onPlay());
popSpeed = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [x0 + w * 0.52 0.82 w * 0.48 0.05], ...
                     'String', compose('%dx real time', speeds), 'FontSize', 11);
chkGrid = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [x0 0.76 w 0.04], ...
                    'String', 'Grid lines', 'Value', 1, 'FontSize', 11, 'BackgroundColor', 'w', 'Callback', @(~, ~) draw());
chkBad = uicontrol(fig, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [x0 0.72 w 0.04], ...
                   'String', 'Points (red: failed)', 'Value', 1, 'FontSize', 11, 'BackgroundColor', 'w', 'Callback', @(~, ~) draw());
info = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.3 w 0.4], 'FontSize', 12, ...
                 'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'FontName', 'FixedWidth');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.04 w 0.24], 'FontSize', 10, ...
          'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'ForegroundColor', [0.35 0.35 0.35], 'String', { ...
          'Green dots: points tracked in this frame.', 'Red crosses: points that failed here and were filled from their neighbors.', ...
          'Blue lines: every 4th grid row and column, deforming with the tissue.', 'Dashed: the tine-tip rectangle at rest.', '', ...
          'Valid points should stay near 100%. A patch of red that persists, or grid lines that kink, marks a tracking problem.'});

if isfield(cfg, 'tine_tips_px')
    z = cfg.tine_tips_px.z; t = cfg.tine_tips_px.theta;
    set(hTips, 'XData', z([1 2 2 1 1]) + 1, 'YData', t([1 1 2 2 1]) + 1);
end
if nargin > 1 && any(strcmp(names, startName)), set(popDef, 'Value', find(strcmp(names, startName))); end
loadDef();
if nargin > 2
    k = min(max(1, round(startFrame)), numel(trk.times)); set(slider, 'Value', k); draw();
end


    function loadDef()
        playing = false; set(btnPlay, 'String', 'Play');
        name = names{get(popDef, 'Value')};
        trk = load(fullfile(testRoot, 'DIC_results', [name '.mat']));
        imgDir = fullfile(testRoot, 'Images', name);
        T = numel(trk.times);
        set(slider, 'Max', max(T, 2), 'Value', 1, 'SliderStep', [1, max(1, round(0.03 * (T - 1)))] / max(T - 1, 1), 'Enable', util.onOff(T > 1));
        k = 1;
        im = imageOps.readGray(fullfile(imgDir, [trk.frames{1} '.jpg']));
        set(hImg, 'CData', im); set(ax, 'CLim', [0 255]);
        axis(ax, 'image');
        if isfield(cfg, 'tine_tips_px')               % zoomed on the tine region, with a margin
            z = cfg.tine_tips_px.z + 1; t = cfg.tine_tips_px.theta + 1;
            m = 0.35 * max(diff(z), diff(t));
            xlim(ax, [max(0.5, z(1) - m), min(size(im, 2) + 0.5, z(2) + m)]);
            ylim(ax, [max(0.5, t(1) - m), min(size(im, 1) + 0.5, t(2) + m)]);
        else
            xlim(ax, [0.5 size(im, 2) + 0.5]); ylim(ax, [0.5 size(im, 1) + 0.5]);
        end
        draw();
    end

    function onSlider()
        k = round(get(slider, 'Value'));
        draw();
    end

    function onPlay()
        playing = ~playing;
        set(btnPlay, 'String', util.ternary(playing, 'Pause', 'Play'));
        if playing && k >= numel(trk.times), k = 1; end
        while playing && isvalid(fig) && k < numel(trk.times)
            t0 = tic;
            sp = speeds(get(popSpeed, 'Value'));
            dtFrame = median(diff(trk.times));
            k = min(k + max(1, round(sp * 0.1 / dtFrame)), numel(trk.times));   % about 10 updates per second
            set(slider, 'Value', k);
            draw();
            pause(max(0.01, max(1, round(sp * 0.1 / dtFrame)) * dtFrame / sp - toc(t0)));
        end
        if isvalid(fig), playing = false; set(btnPlay, 'String', 'Play'); end
    end

    function draw()
        if isempty(trk) || ~isvalid(fig), return; end
        set(hImg, 'CData', imageOps.readGray(fullfile(imgDir, [trk.frames{k} '.jpg'])));
        x = trk.posx(k, :)'; y = trk.posy(k, :)'; v = trk.valid(k, :)';
        if get(chkBad, 'Value')
            set(hGood, 'XData', x(v) + 1, 'YData', y(v) + 1);
            set(hBad, 'XData', x(~v) + 1, 'YData', y(~v) + 1);
        else
            set([hGood hBad], 'XData', NaN, 'YData', NaN);
        end
        if get(chkGrid, 'Value')
            [gx, gy] = trackGridLines(x + 1, y + 1, trk.shape, LINE_EVERY);
            set(hGrid, 'XData', gx, 'YData', gy);
        else
            set(hGrid, 'XData', NaN, 'YData', NaN);
        end
        F = dicCore.affineFit(trk.X0(v, :), [x(v), y(v)]);
        C = F' * F;
        set(info, 'String', { ...
            sprintf('Frame      %d of %d', k, numel(trk.times)), ...
            sprintf('Time       %.1f s', trk.times(k)), '', ...
            sprintf('lambda_z      %.4f', sqrt(C(1, 1))), ...
            sprintf('lambda_theta  %.4f', sqrt(C(2, 2))), '', ...
            sprintf('Valid      %.1f%%', 100 * mean(v)), ...
            sprintf('ZNCC       %.4f (median)', median(trk.zncc(k, :), 'omitnan')), '', ...
            sprintf('Keyframe   %s', util.ternary(any(trk.keyframes == k), 'yes (new reference)', 'no'))});
        drawnow limitrate
    end

    function onClose(~, ~)
        playing = false;
        delete(fig);
    end

    function label(y, s)
        uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 y w 0.03], 'String', s, ...
                  'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
    end
end
