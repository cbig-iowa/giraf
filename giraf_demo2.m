%% GIRAF demo 1: Noisy & undersampled recovery
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
% in Fourier domain in the noiseless setting
%
% Greg Ongie 4/10/2017
%% load data
load('reproduce/data/dataset_SL.mat','x0'); %load Shepp-Logan phantom data
load('reproduce/data/samples_dataset_SL_usf0p65.mat','sampmask'); %random sampling mask
res = size(x0);  %output resolution 
%% define index sets
ind_samples = find(sampmask);
[A,At] = defAAt(ind_samples,res); %undersampling operators
b0 = A(x0);
rng(1);  %fix random seed for reproducible result
sig = 2; %noise level
b = b0 + sig*(randn(size(b0))+1i*randn(size(b0))); %add complex Gaussian white noise
xinit = At(b);
sampSNR = -20*log10(norm(b-b0)/norm(b0));
fprintf('noise added to samples with SNR=%2.2f\n',sampSNR);
%% global problem settings
settings.filter_siz = [25 25];
settings.res = size(x0);
settings.weighting = 'grad'; %gradient weighting--see "get_kspace_weights.m" for other options
settings.exit_tol = 1e-5;  %exit algorithm if relative change in NRMSE between iterates less than this
settings.lambda = 3e7; %regularization parameter (0=equality constrained)
                       %typically in range (1e6-1e8) for noisy data
                       %needs to be tuned for optimal results
settings.p = 0; %Schatten p penalty (0 <= p <= 1)
%% GIRAF  parameters
param.iter = 50; %number of IRLS iterations
param.eps0 = 0; %inital epsilon (0=auto-init) 
param.eta = 1.5; %epsilon decrease factor (typically between 1.1-2);
param.epsmin = 1e-9;
param.ADMM_iter = 200;
param.ADMM_tol = 1e-4;
param.delta = 10; %ADMM conditioning parameter (typically between 10-1000);
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