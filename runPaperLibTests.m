function results = runPaperLibTests()
%RUNPAPERLIBTESTS Run the unit tests of the library (about one minute).
%   results = runPaperLibTests();  table(results)
root = fileparts(mfilename('fullpath'));
addpath(root); setupPaperLib();
results = runtests(fullfile(root,'tests'));
disp(table(results));
fprintf('%d passed, %d failed, %d incomplete.\n', ...
    nnz([results.Passed]),nnz([results.Failed]),nnz([results.Incomplete]));
end
