%% This make.m is to provides a mex interface to the projection problem
%%   in the top-k multiclass SVM model under Mac or Linux.

mex CXXFLAGS='$CXXFLAGS -O3' ./topkprojectionmex/proj_seminewtonmex.cpp
mex CXXFLAGS='$CXXFLAGS -O3' ./topkprojectionmex/proj_knapsackmex.cpp
mex CXXFLAGS='$CXXFLAGS -O3' ./topkprojectionmex/topksvm_dualfvmex.cpp

% *** libsdca is the package downloaded from the website:
% ***            https://github.com/mlapin/libsdca
% *** which is a package provided by Dr. Lapin.

if(exist('matsdca_prox.mexmaci64', 'file') ~= 3 &&...
            exist('matsdca_prox.mexa64', 'file') ~= 3)
        run('./libsdca/src/make')
        copyfile('./libsdca/src/matsdca_*.mex*64', './');
end
