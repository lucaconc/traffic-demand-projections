%% run_parma_overrides.m
% Parma with the corrected sensor-to-edge matching, plus checks and diagnostics.
% Edit the paths if they differ from the defaults of paperConfig, then run.
% All outputs go to <paper>/matlab/parma_override_run (outside the library).

lib = fullfile(fileparts(fileparts(mfilename('fullpath'))));      % this library
clear functions
addpath(lib); setupPaperLib();

cfg = paperConfig();
out = runParmaWithOverrides( ...
    'OverrideCSV', fullfile(cfg.paperRoot,'matlab','parma_edge_overrides.csv'), ...
    'UseParallel', true);
disp(out.decomposition)
