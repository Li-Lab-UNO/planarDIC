function reviewBranches(testRoot, settings)
% Loading branch of every deformation: from the Stretch start to the rake peak, found automatically. With
% settings.review_branches (desktop only), the branches are shown shaded on the rake forces, and wrong
% ones are set by clicking their start and end. Branches set by hand are saved in loading_branches.json
% in the test folder and used by the tracking and the export. Only the analyzed deformations
% (settings.deformations) are reviewed; DIC_results/loading_branches.png shows every tracked one.

P = util.phaseTiming(testRoot);
D = util.dataTable(testRoot);
f = fullfile(testRoot, 'loading_branches.json');
st = P(strcmp({P.Phase}, 'Stretch'));
info = struct('key', {}, 'name', {}, 'a0', {}, 'a1', {}, 'tRec', {}, 'w0', {}, 'w1', {}, 't0', {}, 't1', {});
allDefs = unique([settings.track_deformations, settings.deformations]);
for n = allDefs
    [a0, tRec] = util.stretchPhase(P, n);
    a1 = util.rakePeak(D, a0, tRec, 2.0);
    [t0, t1] = util.loadingBranch(testRoot, n, P, D);
    dur = tRec - a0;
    name = regexprep(st(n).SetName, '\s+save\s+cyc?le$', '', 'ignorecase');   % 'save cycle' is a label only
    info(end + 1) = struct('key', sprintf('Def%02d', n), 'name', name, 'a0', a0, 'a1', a1, 'tRec', tRec, ...
        'w0', max(D.t(1), min(t0, a0) - max(2, 0.15 * dur)), 'w1', min(D.t(end), max(t1, tRec) + max(4, 0.4 * dur)), ...
        't0', t0, 't1', t1); %#ok<AGROW>
end
infoAll = info;
info = infoAll(ismember(allDefs, settings.deformations));

% ---- ask whether the branches are right, and let the wrong ones be clicked (desktop only)
if settings.review_branches && usejava('desktop')
    while true
        fig = plotBranches(D, info, testRoot, '');
        answer = questdlg({'Does each shaded band cover one loading curve,', ...
            'from where the forces start to rise to where they peak?'}, 'runDIC: loading branches', ...
            'Yes, continue', 'No, fix some', 'Yes, continue');
        close(fig);
        if ~strcmp(answer, 'No, fix some'), break; end
        items = arrayfun(@(b) sprintf('%s  %s%s', b.key, b.name, repmat('   (set by hand)', 1, isHand(b))), info, ...
                         'UniformOutput', false);
        [sel, ok] = listdlg('ListString', items, 'SelectionMode', 'multiple', 'Name', 'runDIC', ...
            'PromptString', 'Loading branches to click again:', 'ListSize', [300 160], 'OKString', 'Click them');
        if ~ok || isempty(sel), continue; end
        for i = sel
            tt = clickBranch(D, info(i));
            if isempty(tt), tt = [info(i).a0 info(i).a1]; end          % back to the automatic branch
            info(i).t0 = tt(1); info(i).t1 = tt(2);
        end
    end
end

% ---- save the branches set by hand
byHand = struct();
for b = info(arrayfun(@isHand, info))
    byHand.(b.key) = [b.t0 b.t1];
    util.log('  %s %s: %.1f-%.1f s, set by hand (automatic: %.1f-%.1f s)', b.key, b.name, b.t0, b.t1, b.a0, b.a1);
end
util.writeJSON(f, struct('reviewed', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'by_hand', byHand));
[~, j] = ismember({info.key}, {infoAll.key}); infoAll(j) = info;     % the reviewed ones, as set
plotBranches(D, infoAll, testRoot, fullfile(testRoot, 'DIC_results', 'loading_branches.png'));
util.log('  %d loading branches, %d set by hand -> loading_branches.json, DIC_results/loading_branches.png', ...
         numel(info), numel(fieldnames(byHand)));
end


function h = isHand(b)
h = abs(b.t0 - b.a0) > 0.05 || abs(b.t1 - b.a1) > 0.05;
end


function tt = clickBranch(D, b)
% Ask for the start and the end of one loading branch; [] to go back to the automatic one.
fig = figure('Name', sprintf('runDIC: loading branch of %s', b.key), 'NumberTitle', 'off', 'Color', 'w', ...
             'WindowState', 'maximized');
ax = axes(fig); hold(ax, 'on');
h = drawForces(ax, D, b, [b.a0 b.a1], [0.8 0.8 0.8]);
legend(h, {'F_z (axial)', 'F_\theta (circumferential)', 'automatic branch'}, 'Location', 'northwest');
head = sprintf('%s %s: click where loading starts, then where it ends (the force peak)', b.key, b.name);
t = []; marks = gobjects(0);
while numel(t) < 2
    title(ax, {head, sprintf(['%d of 2 clicked. Right-click or Backspace removes the last one; ' ...
           'A goes back to the automatic branch'], numel(t))}, 'FontSize', 13);
    [xi, ~, button] = ginput(1);
    if isempty(button), continue; end                       % Enter without a click
    if button == 97 || button == 65                          % A: automatic branch
        close(fig); tt = []; return
    elseif button == 3 || button == 8 || button == 127       % right-click, Backspace or Delete: undo
        if ~isempty(t), t(end) = []; delete(marks(end)); marks(end) = []; end
    elseif button == 1
        xi = round(min(max(xi, b.w0), b.w1) * 10) / 10;      % Data.csv rows are 0.1 s apart
        t(end + 1) = xi; %#ok<AGROW>
        marks(end + 1) = xline(ax, xi, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off'); %#ok<AGROW>
    end
    drawnow
end
title(ax, {head, '2 of 2 clicked'}, 'FontSize', 13); drawnow; pause(0.6); close(fig);
tt = sort(t);
end


function h = drawForces(ax, D, b, span, shade)
% Rake forces around one branch, with `span` shaded and the logged Recover start dotted.
m = D.t >= b.w0 & D.t <= b.w1;
h1 = plot(ax, D.t(m), D.X_Force_mN(m), '-', 'Color', [31 119 180] / 255, 'LineWidth', 1.2);
h2 = plot(ax, D.t(m), D.Y_Force_mN(m), '-', 'Color', [255 127 14] / 255, 'LineWidth', 1.2);
yl = ylim(ax); xlim(ax, [b.w0 b.w1]);
hp = patch(ax, span([1 2 2 1]), yl([1 1 2 2]), shade, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
uistack(hp, 'bottom'); ylim(ax, yl);
xline(ax, b.tRec, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
xlabel(ax, 'time (s)'); ylabel(ax, 'force (mN)'); grid(ax, 'on'); box(ax, 'on');
h = [h1 h2 hp];
end


function fig = plotBranches(D, info, testRoot, path)
% Every branch on the rake forces: shown in a window (returned) with an empty path, else saved to path.
show = isempty(path);
if show
    fig = figure('Name', 'runDIC: loading branches', 'NumberTitle', 'off', 'Color', 'w', 'WindowState', 'maximized', ...
                 'MenuBar', 'none', 'ToolBar', 'none');
else
    d = fileparts(path); if ~isfolder(d), mkdir(d); end
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1000] .* [1 1 1 + (numel(info) > 9) * [0.3 0.25]]);
end
nc = ceil(sqrt(numel(info))); nr = ceil(numel(info) / nc);
for i = 1:numel(info)
    b = info(i);
    ax = subplot(nr, nc, i, 'Parent', fig); hold(ax, 'on');
    if isHand(b), c = [0.95 0.55 0.55]; note = ', set by hand'; else, c = [0.55 0.85 0.55]; note = ''; end
    h = drawForces(ax, D, b, [b.t0 b.t1], c);
    title(ax, sprintf('%s %s%s', b.key, b.name, note), 'Interpreter', 'none');
    if i == 1, legend(h, {'F_z', 'F_\theta', 'loading branch'}, 'Location', 'northwest'); end
end
[~, name] = fileparts(testRoot);
sgtitle(fig, {strrep(name, '_', ' '), 'shaded: loading branch exported (green automatic, red set by hand); dotted: logged Recover start'});
if ~show, print(fig, path, '-dpng', '-r90'); close(fig); end
end
