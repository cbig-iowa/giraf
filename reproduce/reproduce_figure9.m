% Reproduce results in Figure 9 of GIRAF paper:
% Comparison of GIRAF versus the LORAKS
% algorithm for recovery of Fourier data
% under a spatial sparisty assumption.
% Uses the "C-LORAKS" matrix lifting.
% For more details on LORAKS see:
% Haldar, Justin P. "Low-Rank Modeling of Local k-Space Neighborhoods (LORAKS) for Constrained MRI." 
% IEEE transactions on medical imaging 33.3 (2014): 668-681.
%
% To perform comparisons with LORAKS, this code
% requires that you have the LORAKS code installed.
% LORAKS is available for download at:
% http://mr.usc.edu/download/loraks/
%
% Note that LORAKS also requires 'lansvd' from 
% the PROPACK library, which is available here:
% http://sun.stanford.edu/~rmunk/PROPACK/
% or as pre-compiled mex files as part of the SVT
% package:
% http://svt.stanfsord.edu/code.html
%
% After downloading and extracting the code for
% LORAKS and lansvd, add the folders to your
% MATLAB path:

% addpath('etc/LORAKS_implementation');
% addpath('etc/lansvd');

% and run the following patch script:

% apply_loraks_patch;

% This generates the file LORAKS_mod_GO.m
% which is a modified version of LORAKS.m
% that computes the error at every iterate.
% Without this file, you will not be able to
% generate the per iteration error plots
% for LORAKS
%%
addpath(genpath(pwd));
GIRAF_ONLY = true; %set to 'false' if comparing with LORAKS
loraks_compare_pwc_SL_random;
loraks_compare_pwc_SL_structured2;