function s = DIC_settings()
% Settings for runDIC. Edit the values below; the defaults suit the CellScale BioTester tests of this
% project and normally need no change.

% ---- what to run
% Steps, in order: 'info' (sample label and notes), 'thickness' (the thickness window: typed, or measured
% from the photos), 'sort' (images into one folder per cycle, with a completeness check), 'tines' (find
% the tine tips), 'scale' (image scale), 'widths' (tine spans), 'track' (DIC), 'artifacts' (find and
% exclude local tracking errors, then track again if any), 'export' (stress-stretch CSVs, plot and the
% AllData file), 'report' (open the tracking report; it reads the existing tracks, so it can be chosen
% alone).
s.steps = {'info', 'thickness', 'sort', 'tines', 'scale', 'widths', 'track', 'artifacts', 'export', 'report'};

% Label added to the sample name in test_config.json (for example 'healthy'); '' for none. The 'info'
% step asks for it (and for notes) in a window instead.
s.label = '';

% Checks shown in a window (desktop only); false skips the window, e.g. for unattended runs.
s.review_sort = true;                  % the completeness of every cycle after sorting

% Tine tips: 'auto' finds them in the image (if that fails for a rake, you are asked to click its tips);
% 'click' asks you to click all 20 (5 per rake) when the automatic ones are wrong. Clicks are kept as clicked.
s.tine_tips = 'auto';
s.review_tines = true;                 % show the tips found and ask whether they are right

% ---- specimen
% Wall thickness H (mm) for the stresses: each test's own, measured from its photos or typed in the
% thickness window and saved as thickness_mm.H in its test_config.json. Without a thickness the stress
% columns are NaN; the forces, stretches and tensions are still written.
s.ask_thickness = true;                % ask at the start for every test without a thickness (desktop only)
% Thickness measured from the photos (the 'thickness' step): positions every thickness_step_mm along the
% strip, edges searched up to thickness_search_mm either side of the centerline, thickness_trim_mm left
% out at each end of a strip.
s.thickness_step_mm = 0.1;
s.thickness_search_mm = 3;
s.thickness_trim_mm = 0.5;

% ---- deformations
s.reference_deformation = 11;          % its first frame is the common reference of the stretches
s.deformations = 11:19;                % analyzed: the CSVs and the tracking report
s.track_deformations = 11:19;         % tracked (fills in the tracking window); 1:19 adds the preconditioning, only for the AllData file's rows
s.calibration_deformations = [1 19];   % used for the image scale
s.loading_only = true;                 % track each deformation only up to the end of its loading branch
s.ask_tracking = true;                 % ask which deformations to track, and loading only or the whole cycle, after the tine tips (desktop only)
% Export: the deformations and the files written by the export step. Formats: 'csv' (one CSV per
% protocol), 'all_csv' (all_protocols.csv), 'png', 'pdf' and 'fig' (the stress-stretch plot), 'xlsx'
% (the AllData file in CellScale's layout). README.txt is always written.
s.export_deformations = 11:19;
s.export_formats = {'csv', 'all_csv', 'png', 'fig', 'xlsx'};
s.ask_export = true;                   % ask which deformations and files to export, once the tracking is done (desktop only)

% ---- DIC (physical sizes, mm)
s.target_px_per_mm = 30;               % images are resampled (down only) to this scale before tracking
s.subset_mm = 1.0;                     % subset size
s.step_mm = 0.2;                       % grid spacing

% ---- tracking report (viewTracking; DIC_results/tracking_report.csv)
% A deformation is Good when its worst frame has at least tracking_good(1) % of the points tracked and a
% median match quality (ZNCC) of at least tracking_good(2); Problem below tracking_problem; Check between.
s.tracking_good = [93 0.95];
s.tracking_problem = [80 0.90];
s.review_tracking = true;              % open the tracking report at the end of a run that tracked (desktop only)

% ---- computing
% MATLAB's matching already uses several CPU cores, so on a computer with up to about 10 cores the
% deformations are tracked one after another (about 13 min for Def11-Def19 of a test on a 10-core
% computer; tracking several at once was slower there). With more cores, several deformations are tracked
% at once, which needs the Parallel Computing Toolbox: each worker tracks one deformation with its share
% of the cores, and needs about 0.6 GB.
s.use_parallel = 'auto';               % 'auto': at once with the toolbox and at least 16 cores; true or false to choose
s.parallel_workers = 'auto';           % deformations at once: 'auto' (half the cores, at most one per deformation) or a number
s.parallel_threads = 'auto';           % threads per worker: 'auto' shares all the processor's threads (hyperthreads too), or a number
s.progress_s = 1;                      % seconds between tracking progress updates (at least 30 with use_parallel)
end
