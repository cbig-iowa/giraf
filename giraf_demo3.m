%% GIRAF demo 3: Effect of filter size
% Solves the convolutional structured matrix recovery program:
%   min_x ||T(x)||_p^p  s.t. Ax = b
% where P is a sampling operator, and
% T(x) is a block Toeplitz matrix built from the data x
% and ||.||_p^p is the Schatten-p quasi-norm, 0 <= p <= 1.
% GIRAF uses an iterative re-weighted least squares approach, 
% that alternates between:
%     (1) updating an annihilating filter/edge mask, and 
%     (2) solving for the data x that is optimally nulled 
%         by the annhilating filter in a least-squares sense
% This demo shows that larger filter sizes gives improved 
% reconstruction quality with only a modest increase in 
% computation time. This is because at each iteration the 
% GIRAF algorithm only requires an eigendecomposition of an 
% NxN gram matrix, where N is the total number of filter coefficients.
%
% Greg Ongie 4/10/2017
%% load data
load('reproduce/data/dataset_BRAIN.mat','x0'); %load brain phantom data
load('reproduce/data/samples_dataset_BRAIN_usf0p50.mat','sampmask'); %random sampling mask
res = size(x0);  %output resolution 
%% define index sets
ind_samples = find(sampmask);
[A,At] = defAAt(ind_samples,res); %undersampling operators
b = A(x0);
xinit = At(b);
%% global problem settings
settings.res = size(x0);
settings.weighting = 'grad'; %gradient weighting--see "get_kspace_weights.m" for other options
settings.exit_tol = 1e-4;  %exit algorithm if relative change in NRMSE between iterates less than this
settings.lambda = 0; %regularization parameter (0=equality constrained)
settings.p = 0; %Schatten p penalty (0 <= p <= 1)
%% GIRAF  parameters
param.iter = 25; %number of IRLS iterations
param.eps0 = 0; %inital epsilon (0=auto-init) 
param.eta = 1.3; %epsilon decrease factor (typically between 1.1-2.0);
param.epsmin = 1e-9;
param.ADMM_iter = 200;
param.ADMM_tol = 1e-4;
param.delta = 10; %ADMM conditioning parameter (typically between 10-1000);
param.overres = settings.res + 2*settings.filter_siz;
%% run GIRAF
settings.filter_siz = [15 15];
tic;
[x,cost] = giraf(xinit,b,A,At,sampmask,param,settings);
time = toc;
SNR = -20*log10(norm(x(:)-x0(:))/norm(x0(:)));

%display results
figure(1); i=1;
subplot(2,3,i); imagesc(abs(ifft2(x)),[0,1]); colorbar; title(sprintf('GIRAF recon\n15x15 filter, SNR=%2.2f\nCPU time=%3.1f',SNR,time)); axis square
subplot(2,3,i+3); imagesc(abs(ifft2(x-x0))); colorbar; title('error'); axis square; drawnow;
%% run GIRAF
settings.filter_siz = [25 25];
tic;
[x,cost] = giraf(xinit,b,A,At,sampmask,param,settings);
time = toc;
SNR = -20*log10(norm(x(:)-x0(:))/norm(x0(:)));

%display results
figure(1); i=2;
subplot(2,3,i); imagesc(abs(ifft2(x)),[0,1]); colorbar; title(sprintf('GIRAF recon\n25x25 filter, SNR=%2.2f\nCPU time=%3.1f',SNR,time)); axis square
subplot(2,3,i+3); imagesc(abs(ifft2(x-x0))); colorbar; title('error'); axis square; drawnow;
%% run GIRAF
settings.filter_siz = [35 35];
tic;
[x,cost] = giraf(xinit,b,A,At,sampmask,param,settings);
time = toc;
SNR = -20*log10(norm(x(:)-x0(:))/norm(x0(:)));

%display results
figure(1); i=3;
subplot(2,3,i); imagesc(abs(ifft2(x)),[0,1]); colorbar; title(sprintf('GIRAF recon\n35x35 filter, SNR=%2.2f\nCPU time=%3.1f',SNR,time)); axis square
subplot(2,3,i+3); imagesc(abs(ifft2(x-x0))); colorbar; title('error'); axis square; drawnow;