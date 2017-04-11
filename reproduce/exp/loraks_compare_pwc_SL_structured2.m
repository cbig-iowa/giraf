% Head-to-head comparison with LORAKS code
load loraks_data.mat;
load samples_loraks_data_structured2.mat;
name = 'loraks_compare_pwc_SL_structured2';
%%
%global settings
settings.R = 4;
settings.filter_siz = [2*settings.R+1,2*settings.R+1];
settings.res = size(x0);
settings.weighting = 'none';
settings.exit_tol = 1e-5;
settings.cost_computations = false;

ind_samples = find(sampmask~=0);
[A,At] = defAAt(ind_samples,settings.res);

%LORAKS settings
param_loraks.lambda = 1e-10; %regularization parameter (set to 0 for equality constraints)
param_loraks.max_iter = 500;  %max number of iterations
param_loraks.tol = 1e-5;
param_loraks.data_scale = scale;
param_loraks.data = 'fulldata.mat';
param_loraks.mask = 'kMask_structured2.mat';

param_loraks.r = 30;

exp_loraks_r30.name = 'LORAKS, r30';
exp_loraks_r30.opt_handle = @opt_loraks_wrapper;
exp_loraks_r30.settings = settings;
exp_loraks_r30.param = param_loraks;

param_loraks.r = 35;

exp_loraks_r35.name = 'LORAKS, r35';
exp_loraks_r35.opt_handle = @opt_loraks_wrapper;
exp_loraks_r35.settings = settings;
exp_loraks_r35.param = param_loraks;

param_loraks.r = 40;

exp_loraks_r40.name = 'LORAKS, r40';
exp_loraks_r40.opt_handle = @opt_loraks_wrapper;
exp_loraks_r40.settings = settings;
exp_loraks_r40.param = param_loraks;

param_loraks.r = 45;

exp_loraks_r45.name = 'loraks_r45';
exp_loraks_r45.opt_handle = @opt_loraks_wrapper;
exp_loraks_r45.settings = settings;
exp_loraks_r45.param = param_loraks;

%GIRAF settings
param_giraf.lambda = 1e-10;%1e-10; %regularization parameter
param_giraf.iter = 100; %number of iterations
param_giraf.eta = 1.1;%1.3; %epsilon decrease factor;    end   
param_giraf.eps = 1e-2;
param_giraf.epsmin = 1e-9;
param_giraf.ADMM_iter = 200;
param_giraf.ADMM_tol = 1e-4;
param_giraf.overres = settings.res + 2*settings.filter_siz;
param_giraf.p = 0;
param_giraf.gam_fac = 1000;
exp_giraf_p0.name = 'GIRAF-0';
exp_giraf_p0.opt_handle = @opt_giraf_flex; %GIRAF modified for circular LORAKS neighborhoods
exp_giraf_p0.settings = settings;
exp_giraf_p0.param = param_giraf;
%% 
if GIRAF_ONLY
    exp = {exp_giraf_p0};
else
    exp = {exp_loraks_r30,exp_loraks_r35,exp_loraks_r40,exp_loraks_r45,exp_giraf_p0}; 
end

results = run_experiments(x0,b,A,At,sampmask,exp);
plot_convergence_results(results);

%show images
if GIRAF_ONLY
    x_giraf  = results{1}.x;
    figure;
    ax1 = subplot(2,1,1);
    imagesc(abs(ifft2(x_giraf)),[0,1]);
    title('GIRAF-0 result');
    colormap(ax1,gray);
    colorbar;
    axis square;   
    ax2 = subplot(2,1,2);
    imagesc(abs(ifft2(x_giraf-x0)),[0,0.1]);
    title('error');
    colormap(ax2,morgenstemning);
    colorbar;
    axis square;    
else
    x_loraks = results{3}.x;
    x_giraf = results{5}.x;
    figure;
    ax1 = subplot(2,2,1);
    imagesc(abs(ifft2(x_giraf)));
    title('LORAKS,r=40 recon');
    colormap(ax1,gray);
    colorbar;
    axis square;   
    ax2 = subplot(2,2,3);
    imagesc(abs(ifft2(x_giraf-x0)));
    title('LORAKS,r=40 error');
    colormap(ax2, morgenstemning);
    colorbar;
    axis square;    
    ax3 = subplot(2,2,2);
    imagesc(abs(ifft2(x_giraf)));
    title('GIRAF-0 recon');
    colormap(ax3,gray);
    colorbar;
    axis square;   
    ax4 = subplot(2,2,4);
    imagesc(abs(ifft2(x_giraf-x0)));
    title('GIRAF-0 error');
    colormap(ax4,morgenstemning);
    colorbar;
    axis square;
end
%save([name,'_results.mat'],'results');