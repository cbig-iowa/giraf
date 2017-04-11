% Reproduce results in Table 1 of GIRAF paper:
% Computation time comparison of GIRAF
% versus competing algorithms (AP, SVT, SVT+UV, IRLS).
% Algorithms are tested on the recovery of 
% piecewise constant images 
% from undersampled Fourier data 
% using the gradient weighted 
% Toeplitz matrix lifting.
% 
% Note: CPU times will vary depending on computer
% but the iteration counts should match those reported
% in the GIRAF paper.

addpath('exp','etc','algs','data');
GIRAF_ONLY = false; %'true'--run GIRAF-0 only
                    %'false'--run all algorithms
                    
% Run experiments to reproduce results in rows 1 & 2                  
row1_PWC1_usf0p50;
row2_PWC1_usf0p33;

% Uncomment lines below to reproduce full table.
% WARNING: may take hours to run all algorithms.
% Alternatively, set "GIRAF_ONLY = true" to see
% results of only the GIRAF algorithm. 
% This will only take 1-5 minutes to run.

% row3_PWC2_usf0p50;
% row4_PWC2_usf0p33;
% row5_SL_usf0p65;
% row6_SL_usf0p50;
% row7_BRAIN_usf0p65;
% row8_BRAIN_usf0p50;
