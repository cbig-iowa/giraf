% Experiment 1: Recovery of synthetic PWC image from 50% uniform random noiseless 
% Fourier samples.
load dataset_SL.mat;
load samples_exp2_SL_usf0p65.mat;
name = 'noisy_SL_usf0p65';
%%
%global settings
settings.filter_siz = [35,35];
settings.r = r;
settings.res = size(x0);
settings.weighting = 'grad';
settings.exit_tol = 1e-5;    %exit algorithm if NMSE below this 
settings.cost_computations = false;

ind_samples = find(sampmask~=0);
[A,At] = defAAt(ind_samples,settings.res);

%AP-PROX settings
param_ap.r = settings.r;
param_ap.lambda = 1e-2; %regularization parameter (set to 0 for equality constraints)
param_ap.iter = 100;  %number of iterations
param_ap.svd_type = 'econ';     %'econ'  = MATLAB svd with 'econ' option (best for large r)
                                %'lansvd'= Lanczos partial svd (best for small r)
                                %'rsvd'  = Random partial svd (best for small r)
exp_ap_prox.name = 'AP-PROX';
exp_ap_prox.opt_handle = @opt_ap_prox;
exp_ap_prox.settings = settings;
exp_ap_prox.param = param_ap;

%SVT+UV settings
param_svt_uv.r = settings.r; %rank cutoff (inner dimension of U,V matricies)
param_svt_uv.lambda = 1e2; %regularization parameter (set to 0 for equality constraints)
param_svt_uv.beta = 10;   %ADMM parameter
param_svt_uv.iter = 50;  %number of ADMM iterations

exp_svt_uv.name = 'SVT+UV';
exp_svt_uv.opt_handle = @opt_svt_uv;
exp_svt_uv.settings = settings;
exp_svt_uv.param = param_svt_uv;

%GIRAF settings
param_giraf.lambda = 3e7; %regularization parameter
param_giraf.iter = 30; %number of iterations
param_giraf.eta = 2;%1.3 %epsilon decrease factor;
param_giraf.eps = 1e-1;
param_giraf.epsmin = 1e-9;
param_giraf.ADMM_iter = 200;
param_giraf.ADMM_tol = 1e-4;
param_giraf.overres = settings.res + 2*settings.filter_siz;

exp_giraf_p0.name = 'GIRAF-0';
exp_giraf_p0.opt_handle = @opt_giraf;
exp_giraf_p0.settings = settings;
exp_giraf_p0.param = param_giraf;
exp_giraf_p0.param.p = 0;
exp_giraf_p0.param.gam_fac = 100;
%% Run experiments
if ~exist('GIRAF_ONLY','var')
    GIRAF_ONLY = true;
end
if GIRAF_ONLY
    exp = {exp_giraf_p0};
else
    exp = {exp_ap_prox,exp_svt_uv,exp_giraf_p0};
end
results = run_experiments(x0,b,A,At,sampmask,exp);
%print_convergence_results(results,1e-4);
plot_convergence_results(results);
%save([name,'_results.mat'],'results');

%show images
if GIRAF_ONLY
    x_giraf  = results{1}.x;
    figure;
    ax = subplot(2,1,1);
    imagesc(abs(ifft2(x_giraf)),[0 1]);
    title('GIRAF-0 result');
    colormap(ax,gray);
    colorbar;
    axis square;   
    
    ax = subplot(2,1,2);
    imagesc(abs(ifft2(x_giraf-x0)),[0 0.1]);
    title('GIRAF-0 error');
    colormap(ax,morgenstemning);
    colorbar;
    axis square;    
else
    x_ap = results{1}.x;
    x_svt = results{2}.x;
    x_giraf = results{3}.x;
    
    figure;
    ax = subplot(2,3,1);
    imagesc(abs(ifft2(x_ap)),[0 1]);
    title('AP-PROX recon');
    colormap(ax,gray);
    colorbar;
    axis square;   
    
    ax = subplot(2,3,4);
    imagesc(abs(ifft2(x_ap-x0)),[0 0.1]);
    title('AP-PROX error');
    colormap(ax, morgenstemning);
    colorbar;
    axis square;
    
    ax = subplot(2,3,2);
    imagesc(abs(ifft2(x_uvt)),[0 1]);
    title('SVT+UV recon');
    colormap(ax,gray);
    colorbar;
    axis square;   
    
    ax = subplot(2,3,5);
    imagesc(abs(ifft2(x_uvt-x0)),[0 0.1]);
    title('SVT+UV error');
    colormap(ax,morgenstemning);
    colorbar;
    axis square;    
    
    ax = subplot(2,3,3);
    imagesc(abs(ifft2(x_giraf)),[0 1]);
    title('GIRAF-0 recon');
    colormap(ax,gray);
    colorbar;
    axis square;   
    
    ax = subplot(2,3,6);
    imagesc(abs(ifft2(x_giraf-x0)), [0 0.1]);
    title('GIRAF-0 error');
    colormap(ax,morgenstemning);
    colorbar;
    axis square;
end