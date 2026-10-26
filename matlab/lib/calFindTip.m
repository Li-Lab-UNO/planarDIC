function [c, tip, contrast] = calFindTip(g, c, t, horizontal, out)
% Refine the centerline c and tip coordinate t of one tine in image g (0-based px).
% Centerline: the darkest line across the rod, 15-45 px behind the tip. Tip: where the brightness along
% the centerline (+-2 px) crosses half-way from the rod's gray level to the tissue's. Returns NaN for the
% tip when no crossing is found; contrast = tissue minus rod gray level.

across = -12:12;
along = t + out * (15:45);
if horizontal
    [X, Y] = meshgrid(along, c + across);
else
    [Y, X] = meshgrid(along, c + across);
end
cp = mean(imageOps.sampleLinear(g, X, Y, 'zero'), 2)';
[~, i] = min(cp);
if i > 1 && i < numel(cp)
    den = cp(i - 1) - 2 * cp(i) + cp(i + 1);
    if den > 0, i = i + 0.5 * (cp(i - 1) - cp(i + 1)) / den; end
end
c = c + across(1) + (i - 1);

% axial profile along the centerline, from rod into tissue
s = -40:40;
a = t - out * s;
w = -2:2;
if horizontal
    [X, Y] = meshgrid(a, c + w);
else
    [Y, X] = meshgrid(a, c + w);
end
prof = mean(imageOps.sampleLinear(g, X, Y, 'zero'), 1);
rod = median(prof(s >= -40 & s <= -15));
tis = median(prof(s >= 12 & s <= 40));
contrast = tis - rod;
mid = 0.5 * (rod + tis);
if tis > rod, above = prof > mid; else, above = prof < mid; end
j = find(above & s > -15, 1);
if isempty(j) || j == 1
    tip = NaN;
    return
end
f = (mid - prof(j - 1)) / (prof(j) - prof(j - 1));
tip = t - out * (s(j - 1) + f);
end
