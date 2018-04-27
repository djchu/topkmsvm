%
% Make script for matsdca
%

clc,clear;

% Library version and the names for the mex files
LIBSDCA_VERSION = '0.3.1';
MEX_PROX = 'matsdca_prox';
MEX_SOLVE = 'matsdca_fit';

% Uncomment the -v for debugging
if strcmp(computer('arch'), 'win64')
    COMPFLAGS = {'$COMPFLAGS', '-std=c++11', ... % '-v', ...
        '-I./sdca', '-I./matlab', '-DBLAS_MATLAB', ...
        sprintf('-DLIBSDCA_VERSION=''"%s"''', LIBSDCA_VERSION), ...
        sprintf('-DMEX_PROX=''"%s"''', MEX_PROX), ...
        sprintf('-DMEX_SOLVE=''"%s"''', MEX_SOLVE), ...
        };
    COMPFLAGS = sprintf('COMPFLAGS="%s"', sprintf('%s ', COMPFLAGS{:}));
else
    CXXFLAGS = {'\$CXXFLAGS', '-std=c++11', '-O3', ...
        '-I./', '-I./matlab', '-DBLAS_MATLAB', ...
        sprintf('-DLIBSDCA_VERSION=''"%s"''', LIBSDCA_VERSION), ...
        sprintf('-DMEX_PROX=''"%s"''', MEX_PROX), ...
        sprintf('-DMEX_SOLVE=''"%s"''', MEX_SOLVE), ...
        };
    CXXFLAGS = sprintf('CXXFLAGS="%s"', sprintf('%s ', CXXFLAGS{:}));
end
old_pwd = pwd;
try
    cd(fileparts(mfilename('fullpath')));
    if strcmp(computer('arch'), 'win64')
        mex(COMPFLAGS,'-largeArrayDims', './matlab/mex_prox.cpp', ...
            '-output', MEX_PROX);
        mex(COMPFLAGS,'-largeArrayDims', './matlab/mex_solve.cpp', ...
            '-output', MEX_SOLVE, '-lmwblas');
    else
        mex(CXXFLAGS,'-largeArrayDims', './matlab/mex_prox.cpp', ...
            '-output', MEX_PROX);
        mex(CXXFLAGS,'-largeArrayDims', './matlab/mex_solve.cpp', ...
            '-output', MEX_SOLVE, '-lmwblas');
    end
catch me
    disp(getReport(me));
    fprintf('Check the configured compiler with `mex -setup`.\n');
    fprintf('If the problem persists, try the CMake build.\n');
end
cd(old_pwd);
