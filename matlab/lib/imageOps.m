classdef imageOps
% Image operations used by the pipeline, in base MATLAB (no toolboxes). Pixel coordinates are
% 0-based (x = column - 1, y = row - 1) throughout, as in the Python pipeline and the json files.

    methods (Static)

        function g = readGray(path)
            % Gray image (double, 0-255) from a CellScale jpg: 0.299 R + 0.587 G + 0.114 B.
            im = imread(path);
            if size(im, 3) == 3
                im = 0.299 * double(im(:, :, 1)) + 0.587 * double(im(:, :, 2)) + 0.114 * double(im(:, :, 3));
                g = round(im);                          % 8-bit gray, as a decoder's gray output
            else
                g = double(im);
            end
        end

        function out = gaussBlur(img, sigma)
            % Separable Gaussian blur; kernel width 8 sigma + 1 (odd), mirrored borders (reflect-101).
            k = round(sigma * 8 + 1);
            if mod(k, 2) == 0, k = k + 1; end
            r = (k - 1) / 2;
            x = -r:r;
            w = exp(-x .^ 2 / (2 * sigma ^ 2));
            w = w / sum(w);
            [h, wd] = size(img);
            P = img(imageOps.reflectIdx(h, r), imageOps.reflectIdx(wd, r));
            out = conv2(w(:), w, P, 'valid');
        end

        function idx = reflectIdx(n, r)
            % 1-based indices of 1..n padded by r on both sides, mirrored without repeating the edge.
            idx = (1 - r):(n + r);
            idx(idx < 1) = 2 - idx(idx < 1);
            idx(idx > n) = 2 * n - idx(idx > n);
            idx = min(max(idx, 1), n);
        end

        function out = resizeArea(img, f)
            % Downsample by factor f (< 1) with area averaging: each output pixel is the mean of the
            % input area it covers. Output size round(f * input size). The weight matrices are cached.
            persistent key Wy Wx
            [h, w] = size(img);
            k = [h w f];
            if ~isequal(key, k)
                Wy = imageOps.areaMatrix(h, round(f * h), 1 / f);
                Wx = imageOps.areaMatrix(w, round(f * w), 1 / f)';
                key = k;
            end
            out = full(Wy * img * Wx);
        end

        function W = areaMatrix(nIn, nOut, s)
            % Sparse nOut x nIn matrix of area weights: output i covers input [i*s, (i+1)*s).
            I = []; J = []; V = [];
            for i = 0:nOut - 1
                a = i * s; b = min((i + 1) * s, nIn);
                for j = floor(a):min(ceil(b) - 1, nIn - 1)
                    ov = min(b, j + 1) - max(a, j);
                    if ov > 1e-12
                        I(end + 1) = i + 1; J(end + 1) = j + 1; V(end + 1) = ov; %#ok<AGROW>
                    end
                end
            end
            W = sparse(I, J, V, nOut, nIn);
            W = spdiags(1 ./ full(sum(W, 2)), 0, nOut, nOut) * W;
        end

        function use = usableMask(raw, level, dilate, resample)
            % 1 where a pixel can be used for matching, 0 on saturated glints (raw gray >= level, grown by
            % `dilate` px), from the raw gray image. After resampling, a pixel is usable only if fully clear.
            sat = double(raw >= level);
            if dilate > 0
                sat = double(conv2(sat, ones(2 * dilate + 1), 'same') > 0);
            end
            use = 1 - sat;
            if resample ~= 1
                use = double(imageOps.resizeArea(use, resample) > 0.999);
            end
        end

        function v = sampleLinear(img, x, y, border)
            % Bilinear samples at 0-based (x, y). border 'replicate' clamps; 'zero' gives 0 outside.
            [h, w] = size(img);
            if nargin < 4, border = 'replicate'; end
            x0 = floor(x); y0 = floor(y); tx = x - x0; ty = y - y0;
            if strcmp(border, 'zero')
                inside = @(xx, yy) xx >= 0 & xx <= w - 1 & yy >= 0 & yy <= h - 1;
                get = @(xx, yy) inside(xx, yy) .* img(min(max(yy, 0), h - 1) + 1 + min(max(xx, 0), w - 1) * h);
            else
                get = @(xx, yy) img(min(max(yy, 0), h - 1) + 1 + min(max(xx, 0), w - 1) * h);
            end
            v = (1 - ty) .* ((1 - tx) .* get(x0, y0) + tx .* get(x0 + 1, y0)) ...
                + ty .* ((1 - tx) .* get(x0, y0 + 1) + tx .* get(x0 + 1, y0 + 1));
        end

        function v = sampleCubic(img, x, y, nanOutside)
            % Bicubic (cubic convolution) samples at 0-based (x, y), borders replicated. With nanOutside,
            % NaN where the 4 x 4 stencil leaves the image (x < 1, y < 1, x > w-3, y > h-3).
            G = griddedInterpolant(img, 'cubic', 'nearest');
            v = G(y + 1, x + 1);
            if nargin > 3 && nanOutside
                [h, w] = size(img);
                v(~(x >= 1 & y >= 1 & x <= w - 3 & y <= h - 3)) = NaN;
            end
        end

        function v = sampleNearest(img, x, y)
            % Nearest-pixel samples at 0-based (x, y); 0 outside the image.
            [h, w] = size(img);
            xi = round(x); yi = round(y);
            inside = xi >= 0 & xi <= w - 1 & yi >= 0 & yi <= h - 1;
            v = zeros(size(x));
            v(inside) = img(yi(inside) + 1 + xi(inside) * h);
        end
    end
end
