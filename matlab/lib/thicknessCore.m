classdef thicknessCore
% Wall thickness of a strip lying on its cut face in a photo: the width of the tissue band, measured
% across a centerline that the user clicks (the seed). Base MATLAB only.
%
% Along the smoothed seed, every settings.thickness_step_mm, the image is read along the local
% perpendicular (up to thickness_search_mm either side). Each pixel gets a tissue score from its color
% (CIE Lab) against the tissue's and the local background's (see crossing): 1 is tissue, 0 background,
% on the blue mat and on the steel ruler alike. The edge on each side is where the score first drops
% below one half, to a fraction of a pixel. The centerline is then moved to the midpoints between the
% edges and the strip measured again, with its normals smoothed over at least 2 wall thicknesses, so a
% rough seed does not widen the chords. As in thickness/measure_thickness.py, positions are not used
% where an edge is missing, the chord is lopsided (one half 3 times the other) or crosses a neighbor's,
% or near the strip's ends; positions far from the strip's median (3 robust SDs: a loose fiber) are
% flagged as outliers. S.reason says why a position is not used.

    methods (Static)

        function S = measureStrip(rgb, seed, pxmm, settings)
            % rgb: the photo (uint8, H x W x 3); seed: n x 2 clicked points (x, y; 1-based px); pxmm:
            % image scale. S has, per position: x, y (center), nx, ny (unit normal), rp, rm (distance to
            % the edge on the + and - side, px), thickness_mm, ok (edges found), outlier.
            R = settings.thickness_search_mm * pxmm;
            [lab, x0, y0] = thicknessCore.labCrop(rgb, seed, R + 1.5 * pxmm);
            C = thicknessCore.curve(seed, settings.thickness_step_mm * pxmm, pxmm);
            for pass = 1:2
                S = thicknessCore.edges(lab, x0, y0, C, R, pxmm);
                if pass == 1 && sum(S.ok) >= 3                      % recenter between the edges
                    m = (S.rp - S.rm) / 2;
                    m(~S.ok) = interp1(find(S.ok), m(S.ok), find(~S.ok), 'linear', 'extrap');
                    P = [S.x + m .* S.nx, S.y + m .* S.ny];
                    % as in measure_thickness.py, the normals are smoothed over at least 2 wall thicknesses
                    D = median(S.rp(S.ok) + S.rm(S.ok));
                    C = thicknessCore.curve(P, settings.thickness_step_mm * pxmm, pxmm, max(pxmm, 2 * D));
                end
            end
            % as in measure_thickness.py: a chord much longer on one side than the other ran out through a
            % ragged edge, and chords that cross a neighbor's (tight bends) measure the bend, not the wall
            S.reason = repmat({''}, numel(S.ok), 1);           % why a position is not used
            S.reason(~S.ok) = {'no edge'};
            lop = S.ok & abs(S.rp - S.rm) ./ (S.rp + S.rm) >= 0.5;
            S.reason(lop) = {'lopsided'}; S.ok = S.ok & ~lop;
            crs = thicknessCore.crossingChords(S) & S.ok;
            S.reason(crs) = {'crosses a neighbor'}; S.ok = S.ok & ~crs;
            S.thickness_mm = (S.rp + S.rm) / pxmm;
            S.thickness_mm(~S.ok) = NaN;
            trim = S.s < settings.thickness_trim_mm * pxmm | S.s > S.s(end) - settings.thickness_trim_mm * pxmm;
            S.reason(trim & S.ok) = {'strip end'}; S.ok(trim) = false;
            S.outlier = thicknessCore.outliers(S.thickness_mm, S.ok);
            S.method = 'auto';
        end

        function S = measureManual(edgeA, edgeB, pxmm, settings)
            % Thickness between two clicked edges (n x 2 each): positions along their middle line, chord
            % along its normal to each edge.
            n = 200;
            A = thicknessCore.resample(edgeA, n); B = thicknessCore.resample(edgeB, n);
            if norm(A(1, :) - B(1, :)) + norm(A(end, :) - B(end, :)) > norm(A(1, :) - B(end, :)) + norm(A(end, :) - B(1, :))
                B = flipud(B);                                      % clicked in opposite directions
            end
            C = thicknessCore.curve((A + B) / 2, settings.thickness_step_mm * pxmm, pxmm);
            k = numel(C.x);
            S = struct('x', C.x, 'y', C.y, 'nx', C.nx, 'ny', C.ny, 's', C.s, 'rp', nan(k, 1), 'rm', nan(k, 1));
            for i = 1:k
                p = [C.x(i) C.y(i)]; d = [C.nx(i) C.ny(i)];
                ta = thicknessCore.hitPolyline(p, d, edgeA); tb = thicknessCore.hitPolyline(p, d, edgeB);
                t = [ta tb];
                if any(t > 0) && any(t < 0), S.rp(i) = min(t(t > 0)); S.rm(i) = -max(t(t < 0)); end
            end
            S.ok = isfinite(S.rp) & isfinite(S.rm);
            S.reason = repmat({''}, numel(S.ok), 1); S.reason(~S.ok) = {'no edge'};
            S.thickness_mm = (S.rp + S.rm) / pxmm;
            S.thickness_mm(~S.ok) = NaN;
            S.outlier = thicknessCore.outliers(S.thickness_mm, S.ok);
            S.method = 'manual';
        end

        function [pxmm, N, ratio] = scaleFromLines(L1, L2, sizeMm)
            % Scale from a line across one square (length L1 px) and a long line (L2 px) between grid
            % lines: N = whole squares in the long line, ratio = L2 / L1 (should be near N).
            ratio = L2 / L1;
            N = max(1, round(ratio));
            pxmm = L2 / (N * sizeMm);
        end

        function out = outliers(t, ok)
            % Positions more than 3 robust SDs (1.4826 x MAD) from the median.
            out = false(size(t));
            v = t(ok);
            if numel(v) < 5, return; end
            med = median(v); mad = 1.4826 * median(abs(v - med));
            out = ok & abs(t - med) > 3 * max(mad, 1e-6);
        end

        % ------------------------------------------------------------------ internals
        function C = curve(P, step, pxmm, smoothPx)
            % Smooth curve through points P (n x 2), sampled every `step` px along its length and smoothed
            % over smoothPx (default 1 mm), with unit normals and arc length s.
            if nargin < 4, smoothPx = pxmm; end
            P = P(:, 1:2);
            keep = [true; any(diff(P), 2)]; P = P(keep, :);
            s = [0; cumsum(hypot(diff(P(:, 1)), diff(P(:, 2))))];
            si = (0:step:s(end))';
            if numel(P(:, 1)) >= 3
                x = interp1(s, P(:, 1), si, 'pchip'); y = interp1(s, P(:, 2), si, 'pchip');
            else
                x = interp1(s, P(:, 1), si); y = interp1(s, P(:, 2), si);
            end
            w = max(3, round(smoothPx / step));
            x = smoothdata(x, 'movmean', w); y = smoothdata(y, 'movmean', w);
            tx = gradient(x); ty = gradient(y); L = hypot(tx, ty);
            C = struct('x', x, 'y', y, 'nx', -ty ./ L, 'ny', tx ./ L, 's', si);
        end

        function S = edges(lab, x0, y0, C, R, pxmm)
            % Edge on each side of every position of curve C (see the class help).
            k = numel(C.x);
            r = (-R:0.5:R);                                           % along the normal (px)
            off = [-0.05 0 0.05] * pxmm;                              % averaged across 0.1 mm along the strip
            prof = zeros(k, numel(r), 3);
            for j = 1:3
                acc = 0;
                for o = off
                    X = C.x + o * (-C.ny) + r .* C.nx - x0 + 1;
                    Y = C.y + o * C.nx + r .* C.ny - y0 + 1;
                    acc = acc + interp2(lab(:, :, j), X, Y, 'linear', NaN);
                end
                prof(:, :, j) = acc / numel(off);
            end
            mid = abs(r) <= 0.15 * pxmm;
            tis = squeeze(median(reshape(prof(:, mid, :), [], 3), 1, 'omitnan'));     % tissue color (1 x 3)
            far = abs(r) >= R - 0.5 * pxmm;
            bgAll = squeeze(median(reshape(prof(:, far, :), [], 3), 1, 'omitnan'));   % the strip's background
            Call = norm(tis - bgAll);
            S = struct('x', C.x, 'y', C.y, 'nx', C.nx, 'ny', C.ny, 's', C.s, 'rp', nan(k, 1), 'rm', nan(k, 1));
            for i = 1:k
                p = squeeze(prof(i, :, :));                           % numel(r) x 3
                for side = [1 -1]
                    sel = side * r >= 0; rr = side * r(sel); pp = p(sel, :);
                    [rr, o] = sort(rr); pp = pp(o, :);
                    % background: at the end of this side, unless that looks like tissue (the other arm of a
                    % tightly curled strip), then the strip's; then refined next to the edge found, likewise
                    bg = median(p(far & side * r > 0, :), 1, 'omitnan');
                    if ~(norm(bg - tis) >= 0.5 * Call), bg = bgAll; end
                    e = thicknessCore.crossing(rr, pp, tis, bg);
                    if isfinite(e)
                        nb = rr >= e + 0.3 * pxmm & rr <= e + 1.0 * pxmm;
                        b2 = median(pp(nb, :), 1, 'omitnan');
                        if sum(nb) >= 3 && norm(b2 - tis) >= 0.5 * Call, e = thicknessCore.crossing(rr, pp, tis, b2); end
                    end
                    if side > 0, S.rp(i) = e; else, S.rm(i) = e; end
                end
            end
            S.ok = isfinite(S.rp) & isfinite(S.rm);
        end

        function e = crossing(rr, pp, tis, bg)
            % First distance where the tissue score drops below 1/2 (NaN if tissue and background hardly
            % differ, the start is not tissue, or there is no drop). As in thickness/measure_thickness.py,
            % a pixel brighter than the background scores its color distance from the background over the
            % tissue's; a darker one scores its warmth (b*) over the tissue's (darker layers, the strip's
            % shaded side), so dark shadows on the mat do not count.
            e = NaN;
            C = norm(tis - bg);
            if any(~isfinite(bg)) || C < 6, return; end              % Lab difference too small to tell
            dE = sqrt(sum((pp - bg) .^ 2, 2));
            warm = (pp(:, 3) - bg(3)) / max(tis(3) - bg(3), 10);
            f = dE / C;
            dark = pp(:, 1) <= bg(1);
            f(dark) = warm(dark);
            if ~(f(1) > 0.5), return; end
            j = find(f < 0.5, 1);
            if isempty(j), return; end
            e = rr(j - 1) + (0.5 - f(j - 1)) / (f(j) - f(j - 1)) * (rr(j) - rr(j - 1));
        end

        function [lab, x0, y0] = labCrop(rgb, P, margin)
            % CIE Lab (D65, single) of the part of rgb around the points P, with `margin` px around them.
            [h, w, ~] = size(rgb);
            x0 = max(1, floor(min(P(:, 1)) - margin)); x1 = min(w, ceil(max(P(:, 1)) + margin));
            y0 = max(1, floor(min(P(:, 2)) - margin)); y1 = min(h, ceil(max(P(:, 2)) + margin));
            lab = thicknessCore.srgb2lab(rgb(y0:y1, x0:x1, :));
        end

        function lab = srgb2lab(rgb)
            % sRGB (uint8) to CIE Lab, D65 white.
            v = single(rgb) / 255;
            lin = v / 12.92;
            hi = v > 0.04045; lin(hi) = ((v(hi) + 0.055) / 1.055) .^ 2.4;
            M = single([0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505]);
            wp = single([0.95047 1 1.08883]);
            [h, w, ~] = size(lin);
            XYZ = reshape(reshape(lin, [], 3) * M', h, w, 3) ./ reshape(wp, 1, 1, 3);
            f = XYZ .^ (1 / 3);
            lo = XYZ <= 0.008856; f(lo) = 7.787 * XYZ(lo) + 16 / 116;
            lab = cat(3, 116 * f(:, :, 2) - 16, 500 * (f(:, :, 1) - f(:, :, 2)), 200 * (f(:, :, 2) - f(:, :, 3)));
        end

        function c = crossingChords(S)
            % Positions whose chord (edge to edge) crosses the chord of another position nearby.
            k = numel(S.x); c = false(k, 1);
            A = [S.x - S.rm .* S.nx, S.y - S.rm .* S.ny]; B = [S.x + S.rp .* S.nx, S.y + S.rp .* S.ny];
            ok = find(S.ok);
            for a = 1:numel(ok)
                i = ok(a);
                reach = (S.rp(i) + S.rm(i)) * 2;                     % neighbors within two thicknesses
                for j = ok(abs(S.s(ok) - S.s(i)) <= reach & ok ~= i)'
                    if thicknessCore.segmentsCross(A(i, :), B(i, :), A(j, :), B(j, :)), c(i) = true; break; end
                end
            end
        end

        function x = segmentsCross(p1, p2, q1, q2)
            o = @(a, b, c) sign((b(1) - a(1)) * (c(2) - a(2)) - (b(2) - a(2)) * (c(1) - a(1)));
            x = o(p1, p2, q1) * o(p1, p2, q2) < 0 && o(q1, q2, p1) * o(q1, q2, p2) < 0;
        end

        function Q = resample(P, n)
            % n points evenly spaced along polyline P.
            s = [0; cumsum(hypot(diff(P(:, 1)), diff(P(:, 2))))];
            [s, u] = unique(s); P = P(u, :);
            Q = [interp1(s, P(:, 1), linspace(0, s(end), n)'), interp1(s, P(:, 2), linspace(0, s(end), n)')];
        end

        function t = hitPolyline(p, d, Pl)
            % Signed distances along direction d from p to every crossing of polyline Pl.
            t = [];
            for k = 1:size(Pl, 1) - 1
                a = Pl(k, :); b = Pl(k + 1, :); e = b - a;
                M = [d(:), -e(:)];
                if abs(det(M)) < 1e-9, continue; end
                q = M \ (a - p)';
                if q(2) >= 0 && q(2) <= 1, t(end + 1) = q(1); end %#ok<AGROW>
            end
        end
    end
end
