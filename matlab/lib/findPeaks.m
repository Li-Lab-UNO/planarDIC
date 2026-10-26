function [pk, prom] = findPeaks(x, minProminence, widthRange, minDistance)
% Peaks of the vector x (1-based indices), as scipy.signal.find_peaks(x, prominence=minProminence,
% width=widthRange, distance=minDistance): local maxima (middle of a plateau), thinned to at least
% minDistance samples apart (higher peaks first), then kept when their prominence is >= minProminence and
% their width at half prominence lies within widthRange = [lo hi]. Also returns the prominences.

x = x(:)';
n = numel(x);
% local maxima, plateaus counted once at their middle
pk = [];
i = 2;
while i < n
    if x(i - 1) < x(i)
        j = i;
        while j < n && x(j + 1) == x(i), j = j + 1; end
        if j < n && x(j + 1) < x(i)
            pk(end + 1) = floor((i + j) / 2); %#ok<AGROW>
            i = j;
        end
    end
    i = i + 1;
end
if isempty(pk), prom = []; return; end

% distance: keep higher peaks first, drop lower ones closer than minDistance
keep = true(size(pk));
[~, order] = sort(x(pk), 'descend');
for o = order
    if ~keep(o), continue; end
    near = abs(pk - pk(o)) < minDistance;
    near(o) = false;
    keep(near & keep) = false;
end
pk = pk(keep);

% prominence and the bases it is measured from
prom = zeros(size(pk)); lb = prom; rb = prom;
for q = 1:numel(pk)
    p = pk(q);
    i = p; lmin = x(p); lb(q) = p;
    while i >= 1 && x(i) <= x(p)
        if x(i) < lmin, lmin = x(i); lb(q) = i; end
        i = i - 1;
    end
    i = p; rmin = x(p); rb(q) = p;
    while i <= n && x(i) <= x(p)
        if x(i) < rmin, rmin = x(i); rb(q) = i; end
        i = i + 1;
    end
    prom(q) = x(p) - max(lmin, rmin);
end
m = prom >= minProminence;
pk = pk(m); prom = prom(m); lb = lb(m); rb = rb(m);

% width at half prominence, interpolated between samples
w = zeros(size(pk));
for q = 1:numel(pk)
    p = pk(q); h = x(p) - 0.5 * prom(q);
    i = p;
    while lb(q) < i && h < x(i), i = i - 1; end
    left = i;
    if x(i) < h, left = left + (h - x(i)) / (x(i + 1) - x(i)); end
    i = p;
    while i < rb(q) && h < x(i), i = i + 1; end
    right = i;
    if x(i) < h, right = right - (h - x(i)) / (x(i - 1) - x(i)); end
    w(q) = right - left;
end
m = w >= widthRange(1) & w <= widthRange(2);
pk = pk(m); prom = prom(m);
end
