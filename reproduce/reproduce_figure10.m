% Reproduce results in Figure 10 of GIRAF paper
% Recovery of noisy, undersampled PWC images
% Note: CPU times will vary depending on computer

addpath('exp','etc','algs','data');
GIRAF_ONLY = true; %'true'--run GIRAF-0 only
                    %'false'--run all algorithms:
                    %AP,SVT,SVT+UV,IRLS-0,GIRAF-0               
noisy_SL_usf0p65;
