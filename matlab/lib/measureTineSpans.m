function widths = measureTineSpans(testRoot, settings)
% Each rake's tine span in the rest frame of the reference deformation, and the specimen widths the
% stresses use: w_theta = mean(left, right) divides the axial force F_z, w_z = mean(top, bottom) divides
% F_theta. A rake's span is the outer-to-outer tine centerline distance x 5/4 (five tine spacings); the
% centerline of each tine is a robust line through the rod's center 15-150 px behind its tip, evaluated
% 15 px behind the tip (where the tines hold the tissue). Writes widths_mm into test_config.json.

cfgFile = fullfile(testRoot, 'test_config.json');
cfg = util.readJSON(cfgFile);
sc = util.readJSON(fullfile(testRoot, 'scale_calibration.json'));
folder = util.defFolder(testRoot, settings.reference_deformation);
[~, paths] = util.findFrames(fullfile(testRoot, 'Images', folder));
g = imageOps.gaussBlur(imageOps.readGray(paths{1}), 1.0);
span = struct();
for s = fieldnames(cfg.calibration_tines_px)'
    R = cfg.calibration_tines_px.(s{1});
    cl = zeros(size(R.tips, 1), 1);
    for k = 1:size(R.tips, 1)
        c0 = R.tips(k, 1); c = c0; t = R.tips(k, 2);
        for it = 1:4                                % tip position (calFindTip's centerline is not used)
            [c2, t2] = calFindTip(g, c, t, R.horizontal, R.out);
            if ~isfinite(t2), break; end
            c = c2; t = t2;
        end
        cl(k) = centerline(g, c0, t, R.horizontal, R.out);
    end
    span.(s{1}) = round((cl(end) - cl(1)) * 5 / 4 / sc.px_per_mm, 2);
end
w = [round((span.left + span.right) / 2, 2), round((span.top + span.bottom) / 2, 2)];
if ~isfield(cfg, 'widths_mm'), cfg.widths_mm = struct(); end
cfg.widths_mm.tine_span = w;
cfg.widths_mm.tine_span_by_rake = span;
util.writeJSON(cfgFile, cfg);
widths = cfg.widths_mm;
util.log('  tine spans: left %.2f, right %.2f, top %.2f, bottom %.2f mm -> w_theta %.2f, w_z %.2f mm', ...
         span.left, span.right, span.top, span.bottom, w(1), w(2));
end


function cAt = centerline(g, c, t, horizontal, out)
% Rod centerline `15` px behind the tip: robust line through the rod centers 15-150 px behind it.
ds = 15:5:150;
across = -35:35;
cs = zeros(size(ds));
for q = 1:numel(ds)
    a = t + out * ds(q) + out * (-2:2);
    if horizontal, [X, Y] = meshgrid(a, c + across); else, [Y, X] = meshgrid(a, c + across); end
    prof = mean(imageOps.sampleLinear(g, X, Y, 'zero'), 2)';
    cs(q) = c + rodCenter(prof, across);
end
A = [ones(numel(ds), 1), ds(:)];
w = ones(numel(ds), 1);
for it = 1:3
    coef = (A .* w) \ (cs(:) .* w);
    r = cs(:) - A * coef;
    sd = 1.4826 * median(abs(r(w > 0))) + 0.3;
    w = double(abs(r) < 3 * sd);
end
cAt = coef(1) + coef(2) * 15;
end


function pos = rodCenter(p, across)
% Center of a dark rod in an across-profile: midpoint of the best dark band between bright flanks.
cs = [0, cumsum(p)];
mn = @(a, b) (cs(b + 1) - cs(a + 1)) / (b - a);     % mean of p(a+1 .. b), 0-based a..b-1
best = -inf; pos = NaN;
for w = 10:26
    for i = 6:numel(p) - w - 7
        score = 0.5 * (mn(i - 6, i) + mn(i + w, i + w + 6)) - mn(i, i + w);
        if score > best, best = score; pos = i + (w - 1) / 2; end
    end
end
pos = across(1) + pos;
end
