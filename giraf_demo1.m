%% GIRAF demo 1: Noiseless recovery
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
% This demonstrates GIRAF for recovery of data undersampled
% in Fourier domain with noiseless samples.
%
% Greg Ongie 4/10/2017
%% load data
load('reproduce/data/dataset_SL.mat','x0'); %load Shepp-Logan phantom data
load('reproduce/data/samples_dataset_SL_usf0p50.mat','sampmask'); %random sampling mask
res = size(x0);  %output resolution 
%% define index sets
ind_samples = find(sampmask);
[A,At] = defAAt(ind_samples,res); %undersampling operators
b = A(x0);
xinit = At(b);
%% global problem settings
settings.filter_siz = [25 25];
settings.res = size(x0);
settings.weighting = 'grad'; %gradient weighting--see "get_kspace_weights.m" for other options
settings.exit_tol = 1e-4;  %exit algorithm if relative change in NRMSE between iterates less than this
settings.lambda = 0; %regularization parameter (0=equality constrained)
settings.p = 0; %Schatten p penalty (0 <= p <= 1)
%% GIRAF parameters
param.iter = 25; %number of IRLS iterations
param.eps0 = 0; %inital epsilon (0=auto-init) 
param.eta = 1.3; %epsilon decrease factor (typically between 1.1-1.5);
param.epsmin = 1e-9;
param.ADMM_iter = 200;
param.ADMM_tol = 1e-4;
param.delta = 10; %ADMM conditioning parameter (typically between 10-100);
param.overres = settings.res + 2*settings.filter_siz;
%% run GIRAF
[x,cost] = giraf(xinit,b,A,At,sampmask,param,settings);
SNR = -20*log10(norm(x(:)-x0(:))/norm(x0(:)));
%% display results
figure(1);
subplot(2,2,1); imagesc(abs(ifft2(x0)),[0,1]); colorbar; title('original'); axis square
subplot(2,2,2); imagesc(abs(ifft2(At(b))),[0,1]); colorbar; title('zero-filled ifft'); axis square
subplot(2,2,3); imagesc(abs(ifft2(x)),[0,1]); colorbar; title(sprintf('recon, SNR=%2.2f',SNR)); axis square
subplot(2,2,4); imagesc(abs(ifft2(x-x0))); colorbar; title('error'); axis square