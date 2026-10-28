function res = measureThickness(testRoot, settings)
% Measure a test's wall thickness from its photos (<sample>/Thickness, <ID>/test <N>/). For each photo
% the user sets the scale (a line across one grid square or ruler tick, then a long line between grid
% lines; the squares in between are counted), frames the strip (two corners; the view stays on it) and
% clicks along its middle; the edges are found automatically (thicknessCore) and can be refined. One
% strip per photo, or more on request. Returns res with H (mm, mean of every position used, all strips
% of all photos pooled), sd, n, strips, photos and a source line, or [] when cancelled. Writes
% <test>/thickness/: stations.csv (every position, used or why not), summary.csv (per strip, with its
% scale and frame) and an overlay per photo.

res = [];
[photoDir, photos] = util.thicknessPhotos(testRoot);
if isempty(photos), errordlg(sprintf('No photos found for %s', testRoot), 'runDIC'); return; end
[~, tn] = fileparts(testRoot); testName = regexp(tn, '^test\d+', 'match', 'once');

DISPLAY_PX = 2000;                               % displayed image size (px) for the whole photo
sq = getpref('runDIC', 'thickness_square', struct('value', '', 'unit', 'in'));
units = {'mm', 'cm', 'in'}; unitMm = [1 10 25.4];

% ---- state
iPhoto = 0; rgb = []; W = 0; H = 0; small = []; ds = 1;
state = ''; pts = zeros(0, 2); L1 = NaN; sizeMm = NaN; pxmm = NaN; nSq = 0; scaleOk = false;
S = []; edgeA = zeros(0, 2); stripNo = 1; msg = ''; strips = {}; photoStrips = {}; frame = [];
cancelled = false;

% ---- window
fig = figure('Name', sprintf('runDIC: thickness, %s', tn), 'NumberTitle', 'off', 'Color', 'w', ...
             'WindowState', 'maximized', 'MenuBar', 'none', 'ToolBar', 'none', 'CloseRequestFcn', @onClose, ...
             'WindowScrollWheelFcn', @onScroll, 'KeyPressFcn', @onKey, 'WindowButtonMotionFcn', @onMove);
ax = axes(fig, 'Units', 'normalized', 'Position', [0.01 0.02 0.69 0.95]);
hImg = image(ax, zeros(2, 2, 3, 'uint8'), 'ButtonDownFcn', @onClick);
set(hImg, 'UserData', @addClick);
axis(ax, 'image', 'off'); hold(ax, 'on');
disableDefaultInteractivity(ax); ax.Toolbar.Visible = 'off';
set(ax, 'ButtonDownFcn', @onClick);      % overlays let clicks through (PickableParts none); the axes catch the rest
hLine = plot(ax, NaN, NaN, '-o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'w', 'LineWidth', 2, 'PickableParts', 'none');   % the clicks
hOk = plot(ax, NaN, NaN, '-', 'Color', [1 0.2 0.2], 'LineWidth', 0.8, 'PickableParts', 'none');
hOut = plot(ax, NaN, NaN, '-', 'Color', [1 0.6 0], 'LineWidth', 1.2, 'PickableParts', 'none');
hEx = plot(ax, NaN, NaN, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, 'PickableParts', 'none');
hEdge = plot(ax, NaN, NaN, '.', 'Color', [1 0.2 0.2], 'MarkerSize', 6, 'PickableParts', 'none');
hDone = plot(ax, NaN, NaN, '-', 'Color', [0.2 0.85 0.3], 'LineWidth', 0.6, 'PickableParts', 'none');
hFrame = plot(ax, NaN, NaN, '--', 'Color', [1 0 0], 'LineWidth', 1.5, 'PickableParts', 'none');
hPrev = plot(ax, NaN, NaN, '--', 'Color', [1 0 0], 'LineWidth', 1.5, 'PickableParts', 'none');   % follows the cursor

x0 = 0.72; w = 0.26;
txtTitle = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.92 w 0.06], 'FontSize', 11, ...
                     'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
txtIns = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.66 w 0.26], 'FontSize', 10, ...
                   'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
% controls as in MATLAB's own dialogs (pixels), placed by placePanel at the window's current size
lblSize = uicontrol(fig, 'Style', 'text', 'Position', [0 0 40 16], 'String', 'Size', ...
                    'FontSize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
edtSize = uicontrol(fig, 'Style', 'edit', 'Position', [0 0 70 22], ...
                    'String', sq.value, 'Callback', @(~, ~) clearError());   % txtErr is made below
popUnit = uicontrol(fig, 'Style', 'popupmenu', 'Position', [0 0 70 22], ...
                    'String', units, 'Value', max(1, find(strcmp(units, sq.unit), 1)));
txtErr = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.575 w 0.035], 'FontSize', 10, ...
                   'ForegroundColor', [0.8 0.1 0.1], 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
btn = gobjects(6, 1);
for kb = 1:6
    btn(kb) = util.button(fig, '', [0 0], [], 'Visible', 'off');   % width set to the label in render
end
txtStat = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.08 w 0.16], 'FontSize', 10, ...
                    'HorizontalAlignment', 'left', 'BackgroundColor', 'w', 'FontName', 'FixedWidth');
btnView = util.button(fig, 'Back to the crop', [0 14], @(~, ~) wholePhoto());   % sized for its longer label
set(btnView, 'String', 'Whole photo');
txtZoom = uicontrol(fig, 'Style', 'text', 'Position', [0 14 260 18], ...
          'String', 'Scroll to zoom where the pointer is.', 'HorizontalAlignment', 'left', 'BackgroundColor', 'w', ...
          'ForegroundColor', [0.4 0.4 0.4]);

set(fig, 'SizeChangedFcn', @(~, ~) placePanel());
placePanel();
nextPhoto();
uiwait(fig);
if cancelled || ~isvalid(fig), if isvalid(fig), delete(fig); end, return; end
delete(fig);
res = finish();


    % ================================================================== photos and states
    function nextPhoto()
        if iPhoto > 0, photoStrips{iPhoto} = strips; end
        iPhoto = iPhoto + 1;
        if iPhoto > numel(photos), uiresume(fig); return; end
        strips = {}; stripNo = 1; S = []; pts = zeros(0, 2); edgeA = zeros(0, 2); msg = ''; frame = [];
        set(txtTitle, 'String', sprintf('Loading %s ...', photos{iPhoto})); drawnow
        rgb = imread(fullfile(photoDir, photos{iPhoto}));
        if size(rgb, 3) == 1, rgb = repmat(rgb, 1, 1, 3); end
        [H, W, ~] = size(rgb);
        ds = max(1, ceil(max(H, W) / DISPLAY_PX));
        small = rgb(1:ds:end, 1:ds:end, :);
        wholePhoto();
        setState('scale1');
    end

    function clearError()
        set(txtErr, 'String', '');
    end

    function placePanel()
        % The panel's pixel-sized controls at the window's current size: on opening (a maximized window
        % gets its size only once shown) and whenever it is resized.
        fp = getpixelposition(fig); X0 = x0 * fp(3); ySize = 0.62 * fp(4);
        set(lblSize, 'Position', [X0 ySize 40 16]);
        set(edtSize, 'Position', [X0 + 44 ySize - 3 70 22]);
        set(popUnit, 'Position', [X0 + 122 ySize - 3 70 22]);
        for kk = 1:numel(btn)
            pb = get(btn(kk), 'Position'); set(btn(kk), 'Position', [X0, 0.565 * fp(4) - kk * 28, pb(3:4)]);
        end
        pv = get(btnView, 'Position'); set(btnView, 'Position', [X0 14 pv(3:4)]);
        set(txtZoom, 'Position', [X0 + pv(3) + 12 14 260 18]);
    end

    function setState(st)
        state = st; pts = zeros(0, 2);
        if any(strcmp(st, {'seed', 'edgeA'})), edgeA = zeros(0, 2); end
        render();
    end

    function render()
        set(txtTitle, 'String', sprintf('%s (photo %d of %d)', photos{iPhoto}, iPhoto, numel(photos)));
        set([lblSize edtSize popUnit], 'Visible', util.onOff(strcmp(state, 'scale1')));
        set(txtErr, 'String', ''); stat = ''; B = {};
        switch state
            case 'scale1'
                ins = {'Step 1 of 5: Scale, 1 grid spacing', '', ...
                       '1. Click on a grid line, then on the next one, straight across one grid spacing (or one ruler tick).', ...
                       '2. Type the size of that spacing and choose its unit, then press Next.', '', ...
                       'Right-click to remove a marker.'};
                if size(pts, 1) == 2, stat = sprintf('one grid spacing: %.1f px', L1); end
                B = {'Next', @scaleNext; 'Skip this photo', @() nextPhoto()};
            case 'scale2'
                ins = {'Step 2 of 5: Scale check', '', ...
                       '1. Click on a grid line near the strip, then on another grid line as far away as the photo allows.', ...
                       '2. Check the count of grid spacings between your clicks (below) and press Confirm.', '', ...
                       'Right-click to remove a marker.'};
                B = {'Back', @() setState('scale1')}; scaleOk = false;
                if size(pts, 1) == 2
                    L2 = norm(diff(pts));
                    [p, n, ratio] = thicknessCore.scaleFromLines(L1, L2, sizeMm);
                    if abs(ratio - n) > 0.25
                        ins = [ins, {'', sprintf(['%.2f grid spacings between your clicks: not a whole number, so a click is ' ...
                                                  'off a grid line (here or in step 1). Draw again.'], ratio)}];
                        B = [{'Draw again', @() setState('scale2')}; B];
                    else
                        nSq = n; pxmm = p; u = get(popUnit, 'Value'); scaleOk = true;
                        ins = [ins, {'', sprintf('%d grid spacings between your clicks = %g %s (%.2f mm): scale %.2f px/mm.', n, ...
                                                 round(n * sizeMm / unitMm(u), 4), units{u}, n * sizeMm, pxmm)}];
                        B = [{'Confirm', @() setState('frame'); 'Draw again', @() setState('scale2')}; B];
                    end
                end
            case 'frame'
                ins = {'Step 3 of 5: Crop', '', ...
                       '1. Click two opposite corners of an ROI around the strip.', ...
                       '2. Press Use this crop. The view then stays on it.', '', ...
                       'Crop whole strips only; leave out pieces cut off at the edge of the photo. If this photo has a second strip of this test, you can crop it after this one.'};
                B = {'No crop (whole photo)', @() useFrame([])};
                if size(pts, 1) == 2
                    B = [{'Use this crop', @() useFrame(pts); 'Draw again', @() setState('frame')}; B];
                end
            case 'seed'
                ins = {'Step 4 of 5: Trace the centerline', '', ...
                       '1. Click along the middle of the strip, from one end to the other.', ...
                       '2. Press Done (or Enter, or double-click the last point). The edges are then found for you.', '', ...
                       'Right-click to remove a marker.'};
                stat = sprintf('%d click(s)', size(pts, 1));
                B = {'Done', @seedDone; 'Skip this photo', @() nextPhoto()};
            case 'check'
                ins = {'Step 5 of 5: Check the thickness measurements', '', ...
                       'Each red line should cross the strip and end on both of its edges.', '', ...
                       'Orange lines were left out automatically (a loose fiber or a cut end). Gray lines are ones you removed.'};
                stat = stripStats();
                B = {'Yes, use this strip', @() useStrip(false); 'Use it, and crop another strip in this photo', @() useStrip(true); ...
                     'Remove a section', @() setState('cut'); ...
                     'Redraw the centerline', @() setState('seed'); 'Click the edges myself', @() setState('edgeA'); ...
                     'Skip this strip', @skipStrip};
            case 'cut'
                ins = {'Remove a section', '', 'Click the two ends of the stretch to leave out, along the strip. Its lines turn gray.'};
                stat = stripStats();
                B = {'Cancel', @() setState('check')};
            case 'edgeA'
                ins = {'Click the edges yourself, 1 of 2', '', ...
                       'Click along one edge of the strip from end to end, then press Done.', '', ...
                       'Right-click to remove a marker.'};
                B = {'Done', @edgeADone; 'Cancel', @() setState('seed')};
            case 'edgeB'
                ins = {'Click the edges yourself, 2 of 2', '', ...
                       'Now click along the other edge from end to end, then press Done.', '', ...
                       'Right-click to remove a marker.'};
                B = {'Done', @edgeBDone; 'Cancel', @() setState('seed')};
        end
        if ~isempty(msg), stat = sprintf('%s\n\n%s', stat, msg); end
        set(txtIns, 'String', ins); set(txtStat, 'String', stat);
        for k = 1:6
            if k <= size(B, 1)
                set(btn(k), 'String', B{k, 1}, 'Callback', @(~, ~) B{k, 2}(), 'Visible', 'on');
                util.fitButton(btn(k));
            else
                set(btn(k), 'Visible', 'off');
            end
        end
        shown = btn(1:size(B, 1));                    % one width, the widest label's, text on the left: a tidy column
        if ~isempty(shown)
            p = cell2mat(get(shown(:), 'Position'));
            for k = 1:numel(shown), set(shown(k), 'Position', [p(k, 1:2) max(p(:, 3)) p(k, 4)], 'HorizontalAlignment', 'left'); end
        end
        drawOverlays();
    end

    % ================================================================== actions
    function scaleNext()
        v = str2double(strrep(strtrim(get(edtSize, 'String')), ',', '.'));
        if size(pts, 1) < 2, set(txtErr, 'String', 'Click across one grid spacing first.'); return; end
        if ~(isfinite(v) && v > 0), set(txtErr, 'String', 'Enter the size of the grid spacing, for example 0.5.'); return; end
        u = get(popUnit, 'Value');
        sizeMm = v * unitMm(u); L1 = norm(diff(pts));
        setpref('runDIC', 'thickness_square', struct('value', get(edtSize, 'String'), 'unit', units{u}));
        setState('scale2');
    end

    function useFrame(P)
        % P: two opposite corners, or [] for no frame.
        if isempty(P)
            frame = [];
        else
            frame = [max(0.5, min(P(:, 1))), max(0.5, min(P(:, 2))), min(W + 0.5, max(P(:, 1))), min(H + 0.5, max(P(:, 2)))];
        end
        setState('seed'); wholePhoto();
    end

    function seedDone()
        if size(pts, 1) < 2, set(txtErr, 'String', 'Click at least twice along the strip.'); return; end
        set(txtStat, 'String', 'Finding the edges ...'); drawnow
        S = thicknessCore.measureStrip(rgb, pts, pxmm, settings);
        S.removed = false(size(S.ok));
        msg = '';
        if sum(used()) < 10
            msg = 'Fewer than 10 positions found: redraw the centerline, or click the edges yourself.';
        end
        state = 'check'; pts = zeros(0, 2); zoomToStrip(); render();
    end

    function zoomToStrip()
        % Show the strip just measured, with 2 mm around it.
        m = 2 * pxmm;
        xs = [min(S.x) - m, max(S.x) + m]; ys = [min(S.y) - m, max(S.y) + m];
        side = max(diff(xs), diff(ys) * W / H);
        xs = mean(xs) + [-0.5 0.5] * side; ys = mean(ys) + [-0.5 0.5] * side * H / W;
        xs = min(max(xs, 0.5), W + 0.5); ys = min(max(ys, 0.5), H + 0.5);
        set(ax, 'XLim', xs, 'YLim', ys); refreshImage();
    end

    function edgeADone()
        if size(pts, 1) < 2, set(txtErr, 'String', 'Click at least twice along the edge.'); return; end
        edgeA = pts; state = 'edgeB'; pts = zeros(0, 2); render();
    end

    function edgeBDone()
        if size(pts, 1) < 2, set(txtErr, 'String', 'Click at least twice along the edge.'); return; end
        S = thicknessCore.measureManual(edgeA, pts, pxmm, settings);
        S.removed = false(size(S.ok));
        msg = ''; edgeA = zeros(0, 2);
        state = 'check'; pts = zeros(0, 2); zoomToStrip(); render();
    end

    function useStrip(another)
        % Keep the strip; then the next photo, or (another) frame a further strip of this photo, with
        % its scale kept.
        if sum(used()) < 10, set(txtErr, 'String', 'A strip needs at least 10 positions.'); return; end
        strips{end + 1} = struct('photo', photos{iPhoto}, 'strip', stripNo, 'S', S, 'pxmm', pxmm, ...
                                 'scale', sprintf('%d squares of %.4g mm', nSq, sizeMm), 'frame', frame);
        if ~another, nextPhoto(); return; end
        msg = sprintf('Strip saved (%.3f mm). Crop the next strip.', mean(S.thickness_mm(used())));
        stripNo = stripNo + 1; S = []; frame = []; wholePhoto(); setState('frame');
    end

    function skipStrip()
        S = []; nextPhoto();
    end

    function u = used()
        u = S.ok & ~S.outlier & ~S.removed;
    end

    function s = stripStats()
        u = used(); t = S.thickness_mm(u);
        s = sprintf('thickness %.3f +- %.3f mm\n%d positions used, %d flagged, %d removed, %d without edges', ...
                    mean(t), std(t), sum(u), sum(S.outlier & ~S.removed), sum(S.removed), sum(~S.ok));
    end

    % ================================================================== mouse, keys, view
    function onClick(~, ~)
        addClick(ax.CurrentPoint(1, 1:2), get(fig, 'SelectionType'));
    end

    function addClick(p, type)
        % A click at image point p; type 'normal' (left) or 'alt' (right). Also the image's UserData,
        % so scripts can click.
        switch type
            case 'alt', undo(); return
            case 'open', primary(); return                  % double-click: its first click placed the point
            case 'normal'
            otherwise, return                                % shift-click
        end
        switch state
            case {'scale1', 'scale2', 'frame'}
                if size(pts, 1) >= 2, pts = zeros(0, 2); end
                pts(end + 1, :) = p;
                if strcmp(state, 'scale1') && size(pts, 1) == 2, L1 = norm(diff(pts)); end
                render();
            case {'seed', 'edgeA', 'edgeB'}
                pts(end + 1, :) = p; render();
            case 'cut'
                pts(end + 1, :) = p;
                if size(pts, 1) == 2
                    sa = nearestS(pts(1, :)); sb = nearestS(pts(2, :));
                    S.removed = S.removed | (S.s >= min(sa, sb) & S.s <= max(sa, sb));
                    state = 'check'; pts = zeros(0, 2);
                end
                render();
        end
    end

    function s = nearestS(p)
        [~, i] = min(hypot(S.x - p(1), S.y - p(2))); s = S.s(i);
    end

    function undo()
        if ~isempty(pts), pts(end, :) = []; render(); end
    end

    function primary()
        % Double-click or Enter: proceed, as the step's first button does, once its clicks are placed.
        switch state
            case 'scale1', scaleNext();
            case 'scale2', if scaleOk, setState('frame'); end
            case 'frame', if size(pts, 1) == 2, useFrame(pts); end
            case 'seed', seedDone();
            case 'edgeA', edgeADone();
            case 'edgeB', edgeBDone();
        end
    end

    function onKey(~, evt)
        switch evt.Key
            case {'return', 'enter'}
                primary();
            case {'backspace', 'delete'}, undo();
        end
    end

    function onMove(~, ~)
        % Crosshair over the photo; preview from the last click to the cursor: the scale line (with the
        % live square count on the long line), the frame, or the next segment of a centerline or edge.
        if ~isvalid(hPrev), return; end
        p = ax.CurrentPoint(1, 1:2); xl = xlim(ax); yl = ylim(ax);
        inside = p(1) >= xl(1) && p(1) <= xl(2) && p(2) >= yl(1) && p(2) <= yl(2);
        set(fig, 'Pointer', util.ternary(inside, 'crosshair', 'arrow'));
        if isempty(pts), return; end
        if ~inside
            set(hPrev, 'XData', NaN, 'YData', NaN); return
        end
        q = pts(end, :);
        switch state
            case {'scale1', 'scale2'}
                if size(pts, 1) ~= 1, return; end
                set(hPrev, 'XData', [q(1) p(1)], 'YData', [q(2) p(2)]);
                if strcmp(state, 'scale2')
                    set(txtStat, 'String', sprintf('%.2f grid spacings', norm(p - q) / L1));
                else
                    set(txtStat, 'String', sprintf('%.1f px', norm(p - q)));
                end
            case 'frame'
                if size(pts, 1) ~= 1, return; end
                set(hPrev, 'XData', [q(1) p(1) p(1) q(1) q(1)], 'YData', [q(2) q(2) p(2) p(2) q(2)]);
            case {'seed', 'edgeA', 'edgeB'}
                set(hPrev, 'XData', [q(1) p(1)], 'YData', [q(2) p(2)]);
            otherwise
                return
        end
        drawnow limitrate
    end

    function onScroll(~, evt)
        p = ax.CurrentPoint(1, 1:2);
        xl = xlim(ax); yl = ylim(ax);
        if p(1) < xl(1) || p(1) > xl(2) || p(2) < yl(1) || p(2) > yl(2), return; end
        f = 1.25 ^ evt.VerticalScrollCount;
        wNew = min(W, diff(xl) * f); hNew = min(H, diff(yl) * f);
        if wNew < 50, return; end
        a = (p(1) - xl(1)) / diff(xl); b = (p(2) - yl(1)) / diff(yl);
        nx = p(1) - a * wNew; ny = p(2) - b * hNew;
        nx = min(max(nx, 0.5), W + 0.5 - wNew); ny = min(max(ny, 0.5), H + 0.5 - hNew);
        set(ax, 'XLim', [nx nx + wNew], 'YLim', [ny ny + hNew]);
        refreshImage();
    end

    function wholePhoto()
        % The frame around the strips when one is set (with 5% around it), else the whole photo.
        if isempty(frame)
            set(ax, 'XLim', [0.5 W + 0.5], 'YLim', [0.5 H + 0.5]);
        else
            c = [mean(frame([1 3])), mean(frame([2 4]))];
            side = 1.05 * max(frame(3) - frame(1), (frame(4) - frame(2)) * W / H);
            xs = c(1) + [-0.5 0.5] * side; ys = c(2) + [-0.5 0.5] * side * H / W;
            set(ax, 'XLim', min(max(xs, 0.5), W + 0.5), 'YLim', min(max(ys, 0.5), H + 0.5));
        end
        set(btnView, 'String', util.ternary(isempty(frame), 'Whole photo', 'Back to the crop'));
        refreshImage();
    end

    function refreshImage()
        % The whole photo subsampled, or, zoomed in, the visible part at up to full resolution.
        xl = xlim(ax); yl = ylim(ax);
        if diff(xl) > 0.6 * W
            set(hImg, 'CData', small, 'XData', [1, 1 + (size(small, 2) - 1) * ds], 'YData', [1, 1 + (size(small, 1) - 1) * ds]);
        else
            cx = round(xl + [-0.5 0.5] * diff(xl)); cy = round(yl + [-0.5 0.5] * diff(yl));
            cx = [max(1, cx(1)), min(W, cx(2))]; cy = [max(1, cy(1)), min(H, cy(2))];
            st = max(1, ceil(diff(cx) / DISPLAY_PX));
            set(hImg, 'CData', rgb(cy(1):st:cy(2), cx(1):st:cx(2), :), 'XData', [cx(1), cx(1) + floor(diff(cx) / st) * st], ...
                'YData', [cy(1), cy(1) + floor(diff(cy) / st) * st]);
        end
        drawnow limitrate
    end

    function drawOverlays()
        set([hLine hOk hOut hEx hEdge], 'XData', NaN, 'YData', NaN);
        set(hPrev, 'XData', NaN, 'YData', NaN);
        if any(strcmp(state, {'scale1', 'scale2', 'frame', 'seed', 'edgeA', 'edgeB', 'cut'})) && ~isempty(pts)
            set(hLine, 'XData', pts(:, 1), 'YData', pts(:, 2));
        end
        if strcmp(state, 'frame') && size(pts, 1) == 2
            set(hLine, 'XData', pts([1 2 2 1 1], 1), 'YData', pts([1 1 2 2 1], 2));
        end
        if ~isempty(frame) && ~strcmp(state, 'frame')
            set(hFrame, 'XData', frame([1 3 3 1 1]), 'YData', frame([2 2 4 4 2]));
        else
            set(hFrame, 'XData', NaN, 'YData', NaN);
        end
        if strcmp(state, 'edgeB') && ~isempty(edgeA)
            set(hEx, 'XData', edgeA(:, 1), 'YData', edgeA(:, 2));
        end
        if any(strcmp(state, {'check', 'cut'})) && ~isempty(S)
            u = used();
            [cxs, cys] = chords(S, u); set(hOk, 'XData', cxs, 'YData', cys);
            [cxs, cys] = chords(S, S.outlier & ~S.removed); set(hOut, 'XData', cxs, 'YData', cys);
            [cxs, cys] = chords(S, S.removed & S.ok); set(hEx, 'XData', cxs, 'YData', cys);
            k = S.ok & ~S.removed;
            set(hEdge, 'XData', [S.x(k) + S.rp(k) .* S.nx(k); S.x(k) - S.rm(k) .* S.nx(k)], ...
                       'YData', [S.y(k) + S.rp(k) .* S.ny(k); S.y(k) - S.rm(k) .* S.ny(k)]);
        end
        dx = []; dy = [];                                     % strips already measured on this photo
        for k = 1:numel(strips)
            [a, b] = chords(strips{k}.S, strips{k}.S.ok & ~strips{k}.S.outlier & ~strips{k}.S.removed);
            dx = [dx; a]; dy = [dy; b]; %#ok<AGROW>
        end
        if isempty(dx), dx = NaN; dy = NaN; end
        set(hDone, 'XData', dx, 'YData', dy);
    end

    function onClose(~, ~)
        if strcmp(questdlg('Stop measuring? Nothing is saved.', 'runDIC', 'Stop', 'Proceed', 'Proceed'), 'Stop')
            cancelled = true; uiresume(fig);
        end
    end

    % ================================================================== results
    function r = finish()
        r = [];
        list = [photoStrips{:}];
        if isempty(list)
            msgbox('No strip was measured, so no thickness was set.', 'runDIC'); return
        end
        t = []; rows = {}; summ = {};
        for k = 1:numel(list)
            A = list{k}; Sk = A.S; u = Sk.ok & ~Sk.outlier & ~Sk.removed;
            t = [t; Sk.thickness_mm(u)]; %#ok<AGROW>
            flag = repmat({''}, numel(Sk.ok), 1);
            flag(~Sk.ok) = Sk.reason(~Sk.ok); flag(Sk.outlier & ~Sk.removed) = {'outlier'}; flag(Sk.removed) = {'removed'};
            rows = [rows; [repmat({A.photo, A.strip, Sk.method, A.pxmm}, numel(Sk.ok), 1), ...
                    num2cell([Sk.x, Sk.y, Sk.thickness_mm, double(u)]), flag]]; %#ok<AGROW>
            summ(end + 1, :) = {A.photo, A.strip, Sk.method, A.pxmm, A.scale, mat2str(round(A.frame)), sum(u), ...
                                mean(Sk.thickness_mm(u)), std(Sk.thickness_mm(u)), median(Sk.thickness_mm(u))}; %#ok<AGROW>
        end
        nPh = numel(unique(cellfun(@(a) a.photo, list, 'UniformOutput', false)));
        r = struct('H', mean(t), 'sd', std(t), 'n', numel(t), 'strips', numel(list), 'photos', nPh);
        r.source = sprintf('measured from %d photo(s), %d strip(s), %d positions', nPh, numel(list), numel(t));
        out = fullfile(testRoot, 'thickness');
        if ~isfolder(out), mkdir(out); end
        writetable(cell2table(rows, 'VariableNames', {'photo', 'strip', 'method', 'px_per_mm', 'x_px', 'y_px', ...
                   'thickness_mm', 'used', 'flag'}), fullfile(out, 'stations.csv'));
        writetable(cell2table(summ, 'VariableNames', {'photo', 'strip', 'method', 'px_per_mm', 'scale', 'frame_px', 'n_used', ...
                   'mean_mm', 'sd_mm', 'median_mm'}), fullfile(out, 'summary.csv'));
        for k = 1:numel(photoStrips)
            if ~isempty(photoStrips{k}), saveOverlay(photoStrips{k}, fullfile(out, [regexprep(photos{k}, '\.\w+$', '') '_overlay.png'])); end
        end
        util.log('  %s: thickness %.3f +- %.3f mm (%s) -> %s', testName, r.H, r.sd, r.source, out);
    end

    function saveOverlay(pl, path)
        % The photo around its strips (3 mm margin), with the chords used.
        im = imread(fullfile(photoDir, pl{1}.photo));
        xy = cell2mat(cellfun(@(q) [q.S.x, q.S.y], pl(:), 'UniformOutput', false));
        fr = pl{1}.frame;
        if ~isempty(fr), xy = [xy; fr([1 2]); fr([3 4])]; end
        m = 3 * pl{1}.pxmm;
        cx = round([max(1, min(xy(:, 1)) - m), min(size(im, 2), max(xy(:, 1)) + m)]);
        cy = round([max(1, min(xy(:, 2)) - m), min(size(im, 1), max(xy(:, 2)) + m)]);
        st = max(1, ceil(max(diff(cx), diff(cy)) / 1600));
        f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 1000]);
        a = axes(f, 'Position', [0 0 1 0.93]);
        image(a, cx, cy, im(cy(1):st:cy(2), cx(1):st:cx(2), :)); axis(a, 'image', 'off'); hold(a, 'on');
        if ~isempty(fr), plot(a, fr([1 3 3 1 1]), fr([2 2 4 4 2]), '--', 'Color', [1 0.85 0.2], 'LineWidth', 1.5); end
        lab = {};
        for k = 1:numel(pl)
            Sk = pl{k}.S; u = Sk.ok & ~Sk.outlier & ~Sk.removed;
            [cxs, cys] = chords(Sk, u); plot(a, cxs, cys, '-', 'Color', [1 0.2 0.2], 'LineWidth', 0.5);
            text(a, Sk.x(1), Sk.y(1), sprintf(' %d', pl{k}.strip), 'Color', 'y', 'FontSize', 14, 'FontWeight', 'bold');
            lab{end + 1} = sprintf('strip %d: %.3f +- %.3f mm (n = %d)', pl{k}.strip, mean(Sk.thickness_mm(u)), ...
                                   std(Sk.thickness_mm(u)), sum(u)); %#ok<AGROW>
        end
        title(a, sprintf('%s, %.2f px/mm: %s', pl{1}.photo, pl{1}.pxmm, strjoin(lab, '; ')), 'Interpreter', 'none');
        print(f, path, '-dpng', '-r90'); close(f);
    end
end


function [xs, ys] = chords(S, sel)
% Chords (edge to edge) at the selected positions, as one NaN-separated line.
k = find(sel);
xs = [S.x(k) - S.rm(k) .* S.nx(k), S.x(k) + S.rp(k) .* S.nx(k), nan(numel(k), 1)]';
ys = [S.y(k) - S.rm(k) .* S.ny(k), S.y(k) + S.rp(k) .* S.ny(k), nan(numel(k), 1)]';
xs = xs(:); ys = ys(:);
if isempty(xs), xs = NaN; ys = NaN; end
end
