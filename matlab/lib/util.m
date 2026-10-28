classdef util
% Small helpers: files, CellScale logs, json, loading branches, thickness, sample label and notes,
% statistics, Command Window output. Base MATLAB only.

    methods (Static)

        function [times, paths] = findFrames(folder)
            % CellScale jpgs in `folder` (not its subfolders), sorted by time (s, from the file name).
            d = dir(fullfile(folder, '*.jpg'));
            times = []; paths = {};
            for k = 1:numel(d)
                t = util.frameTime(d(k).name);
                if ~isnan(t)
                    times(end + 1) = t; paths{end + 1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
                end
            end
            [times, o] = sort(times);
            paths = paths(o);
        end

        function t = frameTime(name)
            % Time (s) in a CellScale image name, e.g. 'ID_..._test1.000760.200.jpg' -> 760.2.
            tok = regexp(name, '\.(\d{6})\.(\d{3})\.jpg$', 'tokens', 'once');
            if isempty(tok), t = NaN; else, t = str2double(tok{1}) + str2double(tok{2}) / 1000; end
        end

        function f = firstFile(folder, pattern)
            % First file in `folder` matching a pattern (e.g. '*Data.csv'); error if none.
            d = dir(fullfile(folder, pattern));
            if isempty(d), error('runDIC:missing', 'No %s in %s', pattern, folder); end
            f = fullfile(d(1).folder, d(1).name);
        end

        function P = phaseTiming(testRoot)
            % Logs/*_PhaseTiming.csv as a struct array with fields Set, SetName, Cycle, Phase, Start.
            f = util.firstFile(fullfile(testRoot, 'Logs'), '*PhaseTiming.csv');
            T = readtable(f, 'VariableNamingRule', 'preserve', 'TextType', 'char', 'Delimiter', ',');
            P = struct('Set', num2cell(T.Set), 'SetName', T.SetName, 'Cycle', num2cell(T.Cycle), ...
                       'Phase', T.Phase, 'Start', num2cell(T.("StartTime(s)")));
        end

        function D = dataTable(testRoot)
            % *Data.csv: one row per time (duplicates dropped), with time rounded to 0.1 s.
            f = util.firstFile(testRoot, '*Data.csv');
            T = readtable(f, 'VariableNamingRule', 'preserve', 'Encoding', 'latin1', 'Delimiter', ',');
            [~, ia] = unique(T.Time_S, 'first');
            T = T(sort(ia), :);
            D = struct('t', round(T.Time_S * 10) / 10, 'X_Size_um', T.X_Size_um, 'Y_Size_um', T.Y_Size_um, ...
                       'X_Force_mN', T.X_Force_mN, 'Y_Force_mN', T.Y_Force_mN);
        end

        function v = atTimes(D, field, t)
            % Values of a Data.csv column at times t (s, rounded to 0.1); NaN where a time is missing.
            [tf, loc] = ismember(round(t(:) * 10) / 10, D.t);
            v = nan(numel(t), 1);
            v(tf) = D.(field)(loc(tf));
        end

        function s = readJSON(path)
            s = jsondecode(fileread(path));
        end

        function writeJSON(path, s)
            fid = fopen(path, 'w');
            fprintf(fid, '%s\n', jsonencode(s, 'PrettyPrint', true));
            fclose(fid);
        end

        function q = percentile(v, p)
            % Percentile p (0-100) of the finite values of v, linear interpolation between order statistics.
            v = sort(v(isfinite(v(:))));
            if isempty(v), q = NaN; return; end
            pos = 1 + (numel(v) - 1) * p / 100;
            lo = floor(pos); hi = ceil(pos);
            q = v(lo) + (pos - lo) * (v(hi) - v(lo));
        end

        function [idx, d] = knn(Q, P, k)
            % Indices (and distances) of the k nearest rows of P for each row of Q (brute force).
            k = min(k, size(P, 1));
            D2 = (Q(:, 1) - P(:, 1)') .^ 2 + (Q(:, 2) - P(:, 2)') .^ 2;
            [d, idx] = mink(D2, k, 2);
            d = sqrt(d);
        end

        function tEnd = rakePeak(D, t0, tRecover, slack)
            % End of a loading branch: the time the rakes peak (largest X/X0 + Y/Y0), searched from the
            % Stretch start t0 to `slack` s past the logged Recover start (the log can lag the motion).
            m = D.t >= t0 - 1e-6 & D.t <= tRecover + slack;
            t = D.t(m); x = D.X_Size_um(m); y = D.Y_Size_um(m);
            [~, i] = max(x / x(1) + y / y(1));
            tEnd = round(t(i) * 10) / 10;
        end

        function [t0, t1, byHand] = loadingBranch(testRoot, n, P, D)
            % Loading branch of deformation n (s): the Stretch start to the rake peak, or the start and end
            % set by hand in the test's loading_branches.json, when it has one.
            [t0, tRec] = util.stretchPhase(P, n);
            t1 = util.rakePeak(D, t0, tRec, 2.0);
            byHand = false;
            f = fullfile(testRoot, 'loading_branches.json');
            key = sprintf('Def%02d', n);
            if isfile(f)
                B = util.readJSON(f);
                if isfield(B, 'by_hand') && isstruct(B.by_hand) && isfield(B.by_hand, key)
                    v = B.by_hand.(key); t0 = v(1); t1 = v(2); byHand = true;
                end
            end
        end

        function [t0, tRecover] = stretchPhase(P, n)
            % Stretch start and Recover start (s) of the n-th Stretch phase (deformation DefNN).
            st = P(strcmp({P.Phase}, 'Stretch'));
            s = st(n);
            rec = P(strcmp({P.Phase}, 'Recover') & [P.Set] == s.Set & [P.Cycle] == s.Cycle);
            t0 = round(s.Start * 10) / 10;
            tRecover = round(rec(1).Start * 10) / 10;
        end

        function folder = defFolder(testRoot, n)
            % Name of the DefNN_... image folder of deformation n.
            d = dir(fullfile(testRoot, 'Images', sprintf('Def%02d_*', n)));
            d = d([d.isdir]);
            if numel(d) ~= 1, error('runDIC:folder', 'Expected one Images/Def%02d_* folder in %s', n, testRoot); end
            folder = d(1).name;
        end

        function [folder, photos] = thicknessPhotos(testRoot)
            % The thickness photo folder of a test, <sample>/Thickness*/test <N> (or test<N>), and its
            % .jpg photos; '' and {} when there is none.
            folder = ''; photos = {};
            [sampleDir, name] = fileparts(testRoot);
            n = regexp(name, '^test(\d+)', 'tokens', 'once');
            if isempty(n), return; end
            for d = dir(fullfile(sampleDir, 'Thickness*'))'
                for cand = {sprintf('test %s', n{1}), sprintf('test%s', n{1})}
                    f = fullfile(d.folder, d.name, cand{1});
                    if ~isfolder(f), continue; end
                    p = dir(f); p = {p(~[p.isdir]).name};
                    p = p(~cellfun(@isempty, regexpi(p, '\.jpe?g$', 'once')));
                    if ~isempty(p), folder = f; photos = sort(p); return; end
                end
            end
        end

        function [H, src] = thickness(testRoot, ~)
            % Thickness (mm) for the stresses: thickness_mm.H of the test's test_config.json (written by
            % the thickness window, or by the Python pipeline), else NaN.
            H = NaN; src = 'not set';
            [t, s] = util.savedThickness(testRoot);
            if ~isempty(t), H = t.H; src = 'test_config.json'; if ~isempty(s), src = [src ', ' s]; end; end
        end

        function [t, src] = savedThickness(testRoot)
            % The thickness_mm entry of the test's test_config.json ([] if none) and its source ('' if none).
            t = []; src = '';
            f = fullfile(testRoot, 'test_config.json');
            if ~isfile(f), return; end
            cfg = util.readJSON(f);
            if isfield(cfg, 'thickness_mm') && isstruct(cfg.thickness_mm) && isfield(cfg.thickness_mm, 'H') ...
                    && isfinite(cfg.thickness_mm.H)
                t = cfg.thickness_mm;
                if isfield(t, 'source'), src = t.source; end
            end
        end

        function s = defList(defs)
            % 'Def11-Def19', 'Def01, Def11-Def15, Def19': runs of consecutive deformations.
            parts = {}; k = 1;
            while k <= numel(defs)
                j = k; while j < numel(defs) && defs(j + 1) == defs(j) + 1, j = j + 1; end
                if j > k, parts{end + 1} = sprintf('Def%02d-Def%02d', defs(k), defs(j)); %#ok<AGROW>
                else, parts{end + 1} = sprintf('Def%02d', defs(k)); end %#ok<AGROW>
                k = j + 1;
            end
            s = strjoin(parts, ', ');
        end

        function h = button(parent, label, pos, callback, varargin)
            % A push button the size of MATLAB's own dialogs' (listdlg): default font, 22 px tall, as wide as its
            % label. pos = [x y] (pixels, lower left), or [x y w] for a set width.
            h = uicontrol(parent, 'Style', 'pushbutton', 'Units', 'pixels', 'String', label, 'Callback', callback, ...
                          'Position', [pos(1:2) 70 22], varargin{:});
            if numel(pos) > 2, set(h, 'Position', [pos(1:3) 22]); else, util.fitButton(h); end
        end

        function hs = buttonRow(parent, labels, callbacks, xRight, y)
            % Buttons (util.button) side by side, 8 px apart, the last one ending at xRight (pixels).
            hs = gobjects(1, numel(labels)); x = xRight;
            for k = numel(labels):-1:1
                hs(k) = util.button(parent, labels{k}, [0 y], callbacks{k});
                p = get(hs(k), 'Position'); x = x - p(3); p(1) = x; set(hs(k), 'Position', p); x = x - 8;
            end
        end

        function centerRow(hs, xCenter)
            % Move a row of buttons (pixels) sideways so that it is centered on xCenter.
            p = cell2mat(get(hs(:), 'Position'));
            dx = xCenter - (min(p(:, 1)) + max(p(:, 1) + p(:, 3))) / 2;
            for k = 1:numel(hs), set(hs(k), 'Position', p(k, :) + [dx 0 0 0]); end
        end

        function fitButton(h)
            % Width of a push button to its label (at least 70 px), height 22 px; the left edge stays.
            u = get(h, 'Units'); set(h, 'Units', 'pixels');
            e = get(h, 'Extent'); p = get(h, 'Position');
            set(h, 'Position', [p(1:2) max(70, ceil(e(3)) + 24) 22]); set(h, 'Units', u);
        end

        function answer(fig, a)
            % A button of a question window (or closing it): keep the answer in its UserData and let the
            % code waiting on it (uiwait) go on.
            set(fig, 'UserData', a);
            uiresume(fig);
        end

        function cancel(where)
            % Stop the whole run (runDIC catches it): the user chose Cancel the run, or closed a question.
            error('runDIC:cancelled', 'cancelled at %s', where);
        end

        function setThickness(testRoot, H, source, sd)
            % Save the thickness H (mm) in the test's test_config.json as thickness_mm: H, H_sd (when given)
            % and source, as the Python pipeline does; other entries (e.g. its workbook value) are kept.
            f = fullfile(testRoot, 'test_config.json');
            if isfile(f), cfg = util.readJSON(f); else, cfg = struct(); end
            t = struct();
            if isfield(cfg, 'thickness_mm') && isstruct(cfg.thickness_mm), t = cfg.thickness_mm; end
            t.H = H;
            if nargin > 3 && isfinite(sd), t.H_sd = sd; elseif isfield(t, 'H_sd'), t = rmfield(t, 'H_sd'); end
            t.source = source;
            cfg.thickness_mm = t;
            util.writeJSON(f, cfg);
        end

        function s = sampleName(testRoot, label)
            % Sample name of test_config.json, e.g. 'N1757 LFPA, test2 (10/06/2026), labeled healthy'.
            [sampleDir, name] = fileparts(testRoot);
            [~, sample] = fileparts(sampleDir);
            tok = regexp(name, '^(test\d+)_(\d\d)(\d\d)(\d{4})', 'tokens', 'once');
            if isempty(tok), s = sprintf('%s, %s', strrep(sample, '_', ' '), name);
            else, s = sprintf('%s, %s (%s/%s/%s)', strrep(sample, '_', ' '), tok{1}, tok{2}, tok{3}, tok{4}); end
            label = strtrim(label);
            if ~isempty(label)
                if ~startsWith(label, 'labeled', 'IgnoreCase', true), label = ['labeled ' label]; end
                s = sprintf('%s, %s', s, label);
            end
        end

        function [label, notes] = sampleInfo(testRoot)
            % Label (e.g. 'healthy', '' if none) and notes (cell of lines) from the test's test_config.json.
            label = ''; notes = {};
            f = fullfile(testRoot, 'test_config.json');
            if ~isfile(f), return; end
            cfg = util.readJSON(f);
            if isfield(cfg, 'sample')
                base = util.sampleName(testRoot, '');
                if startsWith(cfg.sample, base), label = regexprep(cfg.sample(numel(base) + 1:end), '^,\s*', ''); end
                label = regexprep(label, '^labeled\s+', '');
            end
            if isfield(cfg, 'notes') && ~isempty(cfg.notes), notes = cellstr(cfg.notes(:))'; end
        end

        function setSampleInfo(testRoot, label, notes)
            % Write the label (into the sample name) and the notes to the test's test_config.json; the
            % other entries are kept.
            f = fullfile(testRoot, 'test_config.json');
            if isfile(f), cfg = util.readJSON(f); else, cfg = struct(); end
            cfg.sample = util.sampleName(testRoot, label);
            cfg.notes = notes;
            util.writeJSON(f, cfg);
        end

        function s = ternary(c, a, b)
            % a if c is true, else b.
            if c, s = a; else, s = b; end
        end

        function s = onOff(c)
            % 'on' or 'off', for uicontrol properties.
            s = util.ternary(c, 'on', 'off');
        end

        function n = shortName(folder)
            % The last part of a folder path.
            [~, n] = fileparts(folder);
        end

        function log(fmt, varargin)
            % Print a line to the Command Window right away.
            util.progress(sprintf(fmt, varargin{:}));
        end

        function progress(msg, live)
            % Print msg to the Command Window; with live = true it stays on one line that the next call
            % overwrites (a progress line updated in place), otherwise it ends with a new line.
            persistent n                                    % characters of the live line now shown
            if isempty(n), n = 0; end
            fprintf('%s', repmat(sprintf('\b'), 1, n));
            if nargin > 1 && live
                fprintf('%s', msg); n = numel(msg);
            else
                fprintf('%s\n', msg); n = 0;
            end
            drawnow('limitrate');
        end
    end
end
