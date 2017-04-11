function [x,stats] = opt_loraks_wrapper(x0,b,A,At,sampmask,param,settings)
%Wrapper to be used with LORAKS code
%Requires LORAKS.m and lansvd.m to available to the PATH
%See ../reproduce_figure9.m for more information.
%
%addpath('etc/LORAKS_implementation'); %add LORAKS to path
%addpath('etc/lansvd'); %add lansvd to path

fprintf('Starting LORAKS (r = %d, lambda = %1.2e)\n',param.r,param.lambda);

r_C = param.r; %rank cutoff
lambda = param.lambda; %regularization parameter
max_iter = param.max_iter; %max number of iterations
tol = param.tol; %convergence tolerance
scale = param.data_scale;
R = settings.R; %LORAKS neighborhood radius size
res = settings.res;
N1 = res(1);
N2 = res(2);

load(param.data); %load data
load(param.mask); %load sampling mask

% Generate LORAKS operators
fprintf('Generating LORAKS Operators\n\n');
[P_C, Ph_C, ~, ~, ~, ~, cc, ~, ~, sizeC, ~, ~] = generate_LORAKS_operators(N1, N2, R);
lambda_C = lambda/prod(sizeC);
kdat = kData.*kMask;

% C-based reconstruction
[reconstructedKC,~,stats] = LORAKS_mod_GO(kdat, kMask, P_C, Ph_C, [], [], [], [], lambda_C, [], [], r_C, 0, 0, cc, [], [], tol, max_iter,kData,settings);
x = scale*fft2(fftshift(ifft2(ifftshift(reconstructedKC))));

stats.MSE = [(norm(kdat(:)-kData(:))/norm(kData(:)))^2, stats.MSE];
stats.time.iter = [0, stats.time.iter];

fprintf('Done! Total runtime: %6.1f s\n',sum(stats.time.iter));
fprintf('Final cost: %2.3e\n',stats.cost(end));
fprintf('Final MSE: %2.3e\n',stats.MSE(end));
fprintf('\n');
end


