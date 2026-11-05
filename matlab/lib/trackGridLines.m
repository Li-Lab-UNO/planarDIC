function [gx, gy] = trackGridLines(x, y, shape, stepMm)
% Rows and columns of a tracked point grid about every LINE_MM mm (every 4th at a 0.2 mm grid step), as
% one NaN-separated line (for plot). Points are stored row by row (x fastest); shape = [ny nx]; stepMm is
% the track's grid step.
LINE_MM = 0.8;
every = 4;
if isfinite(stepMm), every = max(1, round(LINE_MM / stepMm)); end
X = reshape(x, shape(2), shape(1)); Y = reshape(y, shape(2), shape(1));     % column j = grid row j
r = 1:every:shape(1); c = 1:every:shape(2);
Xr = [X(:, r); nan(1, numel(r))]; Yr = [Y(:, r); nan(1, numel(r))];      % along the rows
Xc = [X(c, :)'; nan(1, numel(c))]; Yc = [Y(c, :)'; nan(1, numel(c))];    % along the columns
gx = [Xr(:); Xc(:)]; gy = [Yr(:); Yc(:)];
end
