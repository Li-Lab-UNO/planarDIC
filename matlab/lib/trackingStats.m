function R = trackingStats(testRoot, settings, quiet)
% Tracking quality of every tracked deformation: valid points per frame and the median match quality
% (ZNCC) per frame, their worst frame, and a status: 'good', 'check' or 'problem' from the limits
% settings.tracking_good and settings.tracking_problem ([valid % ZNCC]). R.analyzed marks the analyzed
% deformations (settings.deformations; the others are preconditioning, Def01-Def10). R.verdict is one
% line for the Command Window and the results' README.txt. Writes DIC_results/tracking_report.csv.

if nargin < 3, quiet = false; end
d = dir(fullfile(testRoot, 'DIC_results', 'Def*.mat'));
d = d(cellfun(@isempty, regexp({d.name}, '_fields\.mat$', 'once')));
num = cellfun(@(f) sscanf(f, 'Def%d'), {d.name});
d = d(ismember(num, [settings.track_deformations, settings.deformations]));
if isempty(d), error('runDIC:missing', 'No tracks in %s: run the track step first', fullfile(testRoot, 'DIC_results')); end
n = numel(d);
R = struct('name', {}, 'frames', {}, 'valid_median', {}, 'valid_worst', {}, 'zncc_median', {}, 'zncc_worst', {}, ...
           'worst_frame', {}, 'worst_fraction', {}, 'keyframes', {}, 'status', {}, 'analyzed', {});
for i = 1:n
    t = load(fullfile(d(i).folder, d(i).name), 'valid', 'zncc', 'keyframes');
    v = 100 * mean(t.valid, 2);                       % valid points per frame (%)
    z = median(t.zncc, 2, 'omitnan');                 % median match quality per frame
    T = numel(v);
    [vw, kv] = min(v); [zw, kz] = min(z);
    sv = grade(vw, settings.tracking_good(1), settings.tracking_problem(1));
    sz = grade(zw, settings.tracking_good(2), settings.tracking_problem(2));
    k = kv; if sz > sv || (sz == sv && vw >= settings.tracking_good(1)), k = kz; end   % frame of the worse measure
    names = {'good', 'check', 'problem'};
    R(i) = struct('name', d(i).name(1:end - 4), 'frames', T, 'valid_median', median(v), 'valid_worst', vw, ...
        'zncc_median', median(z), 'zncc_worst', zw, 'worst_frame', k, 'worst_fraction', (k - 1) / max(T - 1, 1), ...
        'keyframes', numel(t.keyframes), 'status', names{max(sv, sz) + 1}, ...
        'analyzed', ismember(sscanf(d(i).name, 'Def%d'), settings.deformations));
end
writetable(struct2table(R), fullfile(testRoot, 'DIC_results', 'tracking_report.csv'));

a = [R.analyzed];
if any(a), v = verdict(R(a), 'deformations'); else, v = 'no analyzed deformation tracked'; end
if any(~a)                                            % preconditioning, in brief
    p = verdict(R(~a), 'preconditioning cycles');
    v = sprintf('%s; preconditioning: %s', v, regexprep(p, '^all \d+ preconditioning cycles good.*', 'all good'));
end
R(1).verdict = ['Tracking: ' v];
if any(~strcmp({R.status}, 'good')), R(1).verdict = [R(1).verdict '. Show the tracking report (viewTracking) to your supervisor']; end
[R.verdict] = deal(R(1).verdict);
if ~quiet, util.log('  %s -> DIC_results/tracking_report.csv', R(1).verdict); end
end


function v = verdict(R, what)
% 'all 9 deformations good (worst frame ...)' or the ones that are not, with their worst frame.
st = {R.status};
if all(strcmp(st, 'good'))
    v = sprintf('all %d %s good (worst frame %.1f%% valid, ZNCC %.3f)', numel(R), what, min([R.valid_worst]), min([R.zncc_worst]));
else
    parts = {};
    for i = find(~strcmp(st, 'good'))
        parts{end + 1} = sprintf('%s %s (worst frame %.1f%% valid, ZNCC %.3f)', strtok(R(i).name, '_'), ...
                                 upper(R(i).status), R(i).valid_worst, R(i).zncc_worst); %#ok<AGROW>
    end
    v = strjoin(parts, '; ');
end
end


function g = grade(x, good, problem)
% 0 good, 1 check, 2 problem (higher is better for both measures).
if x >= good, g = 0; elseif x >= problem, g = 1; else, g = 2; end
end
