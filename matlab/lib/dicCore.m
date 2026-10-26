classdef dicCore
% Subset DIC: inverse-compositional Gauss-Newton (IC-GN) matching with an affine subset shape and the
% zero-normalized sum of squared differences (ZNSSD), and Lagrangian tracking of a grid of material points
% through an image sequence against keyframes (as the Python pipeline's dic2d.py).
% Parameters per subset: p = [u ux uy v vx vy] (displacement and its gradients). Pixel coordinates are
% 0-based. All subsets are matched at once (vectorized).

    methods (Static)

        function prm = trackParams(radius, maskSaturated, resample)
            % Tracking settings (px in the resampled images).
            prm = struct('radius', radius, 'zncc_min', 0.75, 'median_thr', 3.0, 'median_eps', 0.1, ...
                'update_strain', 0.04, 'update_zncc', 0.93, 'max_gap', 60, 'max_iter', 100, 'tol', 1e-3, ...
                'blur_sigma', 0.8, 'mask_saturated', logical(maskSaturated), 'sat_level', 240, 'sat_dilate', 2, ...
                'exclude_px', zeros(0, 4), 'resample', resample, 'progress_s', 30);
        end

        function [img, use] = loadFrame(path, prm, rects)
            % Image and usable-pixel mask ([] when nothing is masked) for tracking settings prm. rects are
            % excluded rectangles [x0 y0 x1 y1] in this image's (resampled) px.
            raw = imageOps.readGray(path);
            img = raw;
            if prm.resample ~= 1, img = imageOps.resizeArea(raw, prm.resample); end
            if prm.blur_sigma > 0, img = imageOps.gaussBlur(img, prm.blur_sigma); end
            use = [];
            if prm.mask_saturated
                use = imageOps.usableMask(raw, prm.sat_level, prm.sat_dilate, prm.resample);
            end
            if nargin > 2 && ~isempty(rects)
                if isempty(use), use = ones(size(img)); end
                [h, w] = size(img);
                for q = 1:size(rects, 1)
                    r = rects(q, :);
                    use(max(0, floor(r(2))) + 1:min(ceil(r(4)) + 1, h), max(0, floor(r(1))) + 1:min(ceil(r(3)) + 1, w)) = 0;
                end
            end
        end

        % ------------------------------------------------------------------ reference subsets
        function ref = buildReference(centers, radius, img, use)
            % Reference every subset (side 2*radius+1, centered at `centers`, N x 2) to image img.
            d = -radius:radius;
            [DX, DY] = meshgrid(d, d);
            ref.radius = radius;
            ref.dx = DX(:)'; ref.dy = DY(:)';                     % 1 x S
            ref.B3 = [ones(numel(d) ^ 2, 1), ref.dx(:), ref.dy(:)];
            ref.M6 = [ones(numel(d) ^ 2, 1), ref.dx(:), ref.dy(:), ref.dx(:) .^ 2, ref.dx(:) .* ref.dy(:), ref.dy(:) .^ 2];
            ref.corners = [-radius radius -radius radius; -radius -radius radius radius];
            ref.centers = centers;
            ref.hasW = ~isempty(use);
            N = size(centers, 1);
            S = numel(ref.dx);
            ref.f = zeros(N, S); ref.fx = ref.f; ref.fy = ref.f; ref.df = zeros(N, 1);
            ref.Hinv = zeros(6, 6, N); ref.ok = false(N, 1);
            if ref.hasW, ref.w = zeros(N, S); ref.fraw = zeros(N, S); end
            ref = dicCore.updateReference(ref, img, (1:N)', centers, use);
        end

        function ref = updateReference(ref, img, idx, centers, use)
            % Re-reference subsets idx to image img, centered at `centers` (numel(idx) x 2).
            if isempty(idx), return; end
            [gx, gy] = gradient(img);
            X = centers(:, 1) + ref.dx; Y = centers(:, 2) + ref.dy;
            f = imageOps.sampleCubic(img, X, Y, true);
            fx = imageOps.sampleCubic(gx, X, Y, true);
            fy = imageOps.sampleCubic(gy, X, Y, true);
            ok = all(isfinite(f), 2);
            f(~isfinite(f)) = 0; fx(~isfinite(fx)) = 0; fy(~isfinite(fy)) = 0;
            if ref.hasW
                w = double(imageOps.sampleNearest(use, X, Y) > 0.5);
                ok = ok & mean(w, 2) >= 0.5;
                fx = fx .* w; fy = fy .* w;
                fraw = f;
                f = (f - sum(w .* f, 2) ./ max(sum(w, 2), 1)) .* w;
                ref.w(idx, :) = w; ref.fraw(idx, :) = fraw;
            else
                f = f - mean(f, 2);
            end
            df = sqrt(sum(f .^ 2, 2)) + 1e-12;
            % Gauss-Newton Hessian J'J from three products (J = [fx, fx dx, fx dy, fy, fy dx, fy dy])
            Axx = (fx .^ 2) * ref.M6; Axy = (fx .* fy) * ref.M6; Ayy = (fy .^ 2) * ref.M6;
            map = [1 2 3; 2 4 5; 3 5 6];
            n = numel(idx);
            Hinv = zeros(6, 6, n);
            for q = 1:n
                H = [reshape(Axx(q, map), 3, 3), reshape(Axy(q, map), 3, 3); reshape(Axy(q, map), 3, 3)', reshape(Ayy(q, map), 3, 3)];
                if ok(q) && cond(H) < 1e8
                    Hinv(:, :, q) = inv(H);
                else
                    ok(q) = false; Hinv(:, :, q) = eye(6);
                end
            end
            ref.centers(idx, :) = centers;
            ref.f(idx, :) = f; ref.fx(idx, :) = fx; ref.fy(idx, :) = fy; ref.df(idx) = df;
            ref.Hinv(:, :, idx) = Hinv; ref.ok(idx) = ok;
        end

        % ------------------------------------------------------------------ matching
        function [p, zncc, conv] = icgn(ref, img, p0, maxIter, tol)
            % Match every subset of ref into img, starting from p0 (N x 6). Returns p, ZNCC and convergence.
            p = p0; N = size(p, 1);
            conv = false(N, 1);
            active = ref.ok;
            [h, w] = size(img);
            r = ref.radius;
            for it = 1:maxIter
                idx = find(active);
                if isempty(idx), break; end
                pa = p(idx, :);
                [xw, yw] = dicCore.warp(ref, idx, pa);
                bad = dicCore.leavesImage(ref, idx, pa, h, w);
                g = imageOps.sampleCubic(img, xw, yw, false);
                g(bad, :) = 0;
                if ref.hasW
                    wi = ref.w(idx, :);
                    g = (g - sum(wi .* g, 2) ./ max(sum(wi, 2), 1)) .* wi;
                else
                    g = g - mean(g, 2);
                end
                dg = sqrt(sum(g .^ 2, 2)) + 1e-12;
                e = ref.f(idx, :) - (ref.df(idx) ./ dg) .* g;
                b = [(ref.fx(idx, :) .* e) * ref.B3, (ref.fy(idx, :) .* e) * ref.B3];       % n x 6
                dp = -reshape(pagemtimes(ref.Hinv(:, :, idx), reshape(b', 6, 1, [])), 6, [])';
                p(idx, :) = dicCore.compose(pa, dp);
                step = sqrt(dp(:, 1) .^ 2 + dp(:, 4) .^ 2 + r ^ 2 * sum(dp(:, [2 3 5 6]) .^ 2, 2));
                done = step < tol | bad;
                conv(idx(done & ~bad)) = true;
                p(idx(bad), :) = NaN;
                active(idx(done)) = false;
            end
            % correlation coefficient at the final warp
            pp = p; pp(~isfinite(pp)) = 0;
            [xw, yw] = dicCore.warp(ref, (1:N)', pp);
            g = imageOps.sampleCubic(img, xw, yw, true);
            if ref.hasW
                wc = ref.w .* isfinite(g); g(~isfinite(g)) = 0;
                nn = max(sum(wc, 2), 1);
                fc = (ref.fraw - sum(wc .* ref.fraw, 2) ./ nn) .* wc;
                gc = (g - sum(wc .* g, 2) ./ nn) .* wc;
                zncc = sum(fc .* gc, 2) ./ ((sqrt(sum(fc .^ 2, 2)) + 1e-12) .* (sqrt(sum(gc .^ 2, 2)) + 1e-12));
            else
                g = g - mean(g, 2, 'omitnan');
                dg = sqrt(sum(g .^ 2, 2, 'omitnan')) + 1e-12;
                zncc = sum(ref.f .* g, 2, 'omitnan') ./ (ref.df .* dg);
            end
            zncc(~ref.ok | ~all(isfinite(p), 2)) = NaN;
        end

        function [xw, yw] = warp(ref, idx, p)
            xw = ref.centers(idx, 1) + p(:, 1) + (1 + p(:, 2)) .* ref.dx + p(:, 3) .* ref.dy;
            yw = ref.centers(idx, 2) + p(:, 4) + p(:, 5) .* ref.dx + (1 + p(:, 6)) .* ref.dy;
        end

        function bad = leavesImage(ref, idx, p, h, w)
            % True for subsets whose bicubic stencil leaves the image (the warp is affine: its corners decide).
            cdx = ref.corners(1, :); cdy = ref.corners(2, :);
            x = ref.centers(idx, 1) + p(:, 1) + (1 + p(:, 2)) .* cdx + p(:, 3) .* cdy;
            y = ref.centers(idx, 2) + p(:, 4) + p(:, 5) .* cdx + (1 + p(:, 6)) .* cdy;
            bad = ~all(x >= 1 & y >= 1 & x <= w - 3 & y <= h - 3, 2);
        end

        function pn = compose(p, dp)
            % Inverse-compositional update: W(p) o W(dp)^-1, as parameters.
            a1 = 1 + p(:, 2); b1 = p(:, 3); c1 = p(:, 1); d1 = p(:, 5); e1 = 1 + p(:, 6); f1 = p(:, 4);
            a = 1 + dp(:, 2); b = dp(:, 3); c = dp(:, 1); d = dp(:, 5); e = 1 + dp(:, 6); f = dp(:, 4);
            det = a .* e - b .* d;
            ia = e ./ det; ib = -b ./ det; ic = (b .* f - c .* e) ./ det;
            id = -d ./ det; ie = a ./ det; if_ = (c .* d - a .* f) ./ det;
            pn = [a1 .* ic + b1 .* if_ + c1, a1 .* ia + b1 .* id - 1, a1 .* ib + b1 .* ie, ...
                  d1 .* ic + e1 .* if_ + f1, d1 .* ia + e1 .* id, d1 .* ib + e1 .* ie - 1];
        end

        % ------------------------------------------------------------------ neighbors
        function p = fillFromNeighbors(pts, p, good, k)
            % Replace p(~good, :) by the inverse-distance-weighted mean of the k nearest good points.
            if nargin < 4, k = 8; end
            if all(good) || sum(good) < 3, return; end
            [j, d] = util.knn(pts(~good, :), pts(good, :), k);
            wt = 1 ./ max(d, 1e-6);
            pg = p(good, :);
            out = zeros(sum(~good), size(p, 2));
            for c = 1:size(p, 2)
                v = pg(:, c);
                out(:, c) = sum(reshape(v(j), size(j)) .* wt, 2) ./ sum(wt, 2);
            end
            p(~good, :) = out;
        end

        function pass = medianTest(du, nb, good, epsv, thr)
            % Normalized median test (Westerweel and Scarano 2005) of displacements du (N x 2) against the
            % good ones of the neighbors nb (N x K). Points with fewer than 3 good neighbors pass.
            G = good(nb);
            r = zeros(size(du));
            for c = 1:2
                v = du(:, c);
                U = reshape(v(nb), size(nb)); U(~G) = NaN;
                med = median(U, 2, 'omitnan');
                res = median(abs(U - med), 2, 'omitnan');
                r(:, c) = abs(v - med) ./ (res + epsv);
            end
            pass = sum(G, 2) < 3 | max(r, [], 2) <= thr;
        end

        function rects = shiftRects(rects, X0, x, valid)
            % Move rectangles [x0 y0 x1 y1], given at reference positions X0, with the tissue (positions x).
            if isempty(rects), return; end
            c = [(rects(:, 1) + rects(:, 3)) / 2, (rects(:, 2) + rects(:, 4)) / 2];
            [j, d] = util.knn(c, X0(valid, :), 8);
            wt = 1 ./ max(d, 1e-6);
            u = x(valid, :) - X0(valid, :);
            ux = u(:, 1); uy = u(:, 2);
            du = [sum(reshape(ux(j), size(j)) .* wt, 2), sum(reshape(uy(j), size(j)) .* wt, 2)] ./ sum(wt, 2);
            rects = rects + [du du];
        end

        % ------------------------------------------------------------------ sequence tracking
        function R = trackSequence(paths, X0, prm, times, label, report)
            % Track material points X0 (N x 2) through image `paths`. A point is valid in a frame when IC-GN
            % converges, ZNCC >= zncc_min and its frame-to-frame displacement passes the median test.
            % Failed points get a position interpolated from valid neighbors (flagged invalid) and keep their
            % last good reference. Keyframes are renewed (for the valid points) when the displacement gradient
            % since the keyframe exceeds update_strain, the median ZNCC drops below update_zncc, or after
            % max_gap frames.
            T = numel(paths); N = size(X0, 1);
            posx = nan(T, N); posy = posx; zn = posx; valid = false(T, N);
            rects0 = prm.exclude_px;
            [img, use] = dicCore.loadFrame(paths{1}, prm, rects0);
            posx(1, :) = X0(:, 1); posy(1, :) = X0(:, 2); zn(1, :) = 1;
            ref = dicCore.buildReference(X0, prm.radius, img, use);
            valid(1, :) = ref.ok;
            nb = util.knn(X0, X0, 9); nb = nb(:, 2:end);            % 8 nearest grid neighbors
            kRef = 1; keyframes = 1;
            p = zeros(N, 6); pPrev = p;
            t0 = tic; tLast = 0;
            for k = 2:T
                [img, use] = dicCore.loadFrame(paths{k}, prm);
                rdt = 1;                                            % constant-velocity predictor, scaled by
                if k >= 3                                           % the actual time step (dropped frames)
                    rdt = (times(k) - times(k - 1)) / max(times(k - 1) - times(k - 2), 1e-9);
                end
                pInit = p + (p - pPrev) * rdt;
                [pNew, c, cv] = dicCore.icgn(ref, img, pInit, prm.max_iter, prm.tol);
                good = cv & c >= prm.zncc_min & ref.ok & all(isfinite(pNew), 2);
                pNew(~isfinite(pNew)) = p(~isfinite(pNew));
                x = ref.centers + pNew(:, [1 4]);
                mt = dicCore.medianTest(x - [posx(k - 1, :)', posy(k - 1, :)'], nb, good, prm.median_eps, prm.median_thr);
                good = good & mt;
                bad = ~good;
                if any(bad) && sum(good) >= 3
                    u = dicCore.fillFromNeighbors(X0, x - X0, good);
                    x(bad, :) = X0(bad, :) + u(bad, :);
                    pNew(bad, 1) = x(bad, 1) - ref.centers(bad, 1);
                    pNew(bad, 4) = x(bad, 2) - ref.centers(bad, 2);
                    pNew(bad, [2 3 5 6]) = p(bad, [2 3 5 6]);
                end
                posx(k, :) = x(:, 1); posy(k, :) = x(:, 2); zn(k, :) = c; valid(k, :) = good;
                pPrev = p; p = pNew;
                if any(good), grad = max(abs(pNew(good, [2 3 5 6])), [], 2); else, grad = 0; end
                need = median(grad) > prm.update_strain || median(c, 'omitnan') < prm.update_zncc || k - kRef >= prm.max_gap;
                if need && k < T
                    idx = find(good);
                    vel = x(idx, :) - [posx(k - 1, idx)', posy(k - 1, idx)'];
                    if ~isempty(rects0)       % excluded rectangles follow the tissue into the new keyframe
                        [~, use] = dicCore.loadFrame(paths{k}, prm, dicCore.shiftRects(rects0, X0, x, good));
                    end
                    ref = dicCore.updateReference(ref, img, idx, x(idx, :), use);
                    p(idx, :) = 0; pPrev(idx, :) = 0;
                    pPrev(idx, 1) = -vel(:, 1); pPrev(idx, 4) = -vel(:, 2);
                    kRef = k; keyframes(end + 1) = k; %#ok<AGROW>
                end
                el = toc(t0);
                if el - tLast >= prm.progress_s || k == T
                    report(sprintf('    %s  %d/%d frames  %3.0f%%  %s elapsed, %s left  (valid %.1f%%, ZNCC %.4f)', label, k, T, ...
                                   100 * k / T, hms(el), hms(el / (k - 1) * (T - k)), 100 * mean(good), median(c, 'omitnan')), k < T);
                    tLast = el;
                end
            end
            R = struct('posx', posx, 'posy', posy, 'zncc', zn, 'valid', valid, 'keyframes', keyframes);
        end

        % ------------------------------------------------------------------ kinematics
        function F = affineFit(X, x)
            % Least-squares homogeneous deformation gradient F mapping reference points X to x (both n x 2).
            F = ((X - mean(X, 1)) \ (x - mean(x, 1)))';
        end

        function F = bilinearCenterF(P0, P)
            % F at the center of a 4-node bilinear element; nodes (4 x 2) counter-clockwise.
            dN = 0.25 * [-1 -1; 1 -1; 1 1; -1 1];
            F = (P' * dN) / (P0' * dN);
        end
    end
end


function s = hms(t)
t = round(t);
if t >= 3600, s = sprintf('%dh%02dm', floor(t / 3600), floor(mod(t, 3600) / 60));
else, s = sprintf('%dm%02ds', floor(t / 60), mod(t, 60)); end
end
