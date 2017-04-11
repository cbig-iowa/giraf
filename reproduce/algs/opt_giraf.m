function [x,stats] = opt_giraf(x0,b,A,At,sampmask,param,settings)
%Solves min_x 0.5*||Ax-b||_2^2 + lambda*||T(x)||_p^p
%where A is the measurement operator, and
%T(x) is the block Toeplitz matrix built from k-space data x
%and ||.||_p is the Shatten p-norm, with 0 <= p <= 1.
%This code is meant for benchmarking only
fprintf('Starting GIRAF (p=%d, lambda = %1.2e)\n',param.p,param.lambda);

lambda = param.lambda; %regularization parameter
iter = param.iter; %number of iterations
eta = param.eta; %epsilon decrease factor;
eps = param.eps;
epsmin = param.epsmin;
p = param.p;
q = 1-(p/2);      
ADMM_iter = param.ADMM_iter;
ADMM_tol = param.ADMM_tol;

if(isfield(param,'gam_fac'))
    gam_fac = param.gam_fac;
else
    gam_fac = 100;
end

filter_siz = settings.filter_siz; 
res = settings.res;
filter_siz2 = 2*filter_siz - [1,1]; %squared filter dimensions
if(isfield(param,'overres'))
    overres = param.overres;
else
    overres = res + 2*filter_siz;       %oversampled reconstruction grid
end
k = get_kspace_inds(overres);
ind_full = get_lowpass_inds(k,res);
ind_filter = get_lowpass_inds(k,filter_siz);
ind_filter2 = get_lowpass_inds(k,filter_siz2);

dz = get_kspace_weights(settings.weighting,k,overres,settings);
 
M = @(z) repmat(z,[1,1,size(dz,3)]).*dz;
Mt = @(Z) sum(Z.*conj(dz),3);
MtMmask = Mt(M(ones(overres))); 

sampmask_pad = zeros(overres);
sampmask_pad(ind_full) = sampmask;

Atb = At(b);
Atb_pad = zeros(overres);
Atb_pad(ind_full) = Atb;

x = Atb;
x_pad = zeros(overres);
x_pad(ind_full) = x;

stats.time.svd_step = zeros(1,iter+1);
stats.time.iter = zeros(1,iter+1);
stats.MSE = zeros(1,iter+1);
stats.cost = zeros(1,iter);

stats.time.svd_step(1) = 0;
stats.time.iter(1) = 0;
stats.MSE(1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;
for i=1:iter
    tic;
    
    %step 1: Compute sos annihilating polynomial
    gradx = M(x_pad);
    G = build_gram_matrix(gradx,filter_siz,filter_siz2,ind_filter2);
    tic;
    [U,S] = eig(G+eps*eye(size(G)));
    stats.time.svd_step(i+1) = toc;
    s = diag(S);
    
    if(settings.cost_computations)    
        if p == 0
            shatten = 0.5*sum(log(s-eps));
        else
            shatten = (1/p)*sum((s-eps).^(p/2));
        end        
        diff = A(x)-b;
        if(lambda == 0) %equality constrained
            stats.cost(i) = shatten;
        else
            stats.cost(i) = 0.5*norm(diff(:)).^2 + lambda*shatten;          
        end  
    end
    mu = build_sos_poly(U,s,q,overres,filter_siz,filter_siz2,ind_filter,ind_filter2);

    %step 2: ADMM solution of least squares problem
    gam = max(mu(:))/gam_fac;  %set ADMM parameter some fraction of sos-mask max value
    [x_pad,~] = run_ADMM_WL2(x_pad,mu,M,Mt,sampmask_pad,MtMmask,Atb_pad,gam,lambda,ADMM_iter,ADMM_tol);
    %[x_pad,flag,relres,iter,resvec] = run_LSQR_WL2(b,mu,overres,k,dz,sampmask_pad,LSQR_iter,lambda);
    x_pad = reshape(x_pad,overres);
    x = reshape(x_pad(ind_full),res);
    
    %update epsilon
    eps = max(eps/eta,epsmin);
    
    stats.time.iter(i+1) = toc; 

    %cost/MSE computations
    stats.MSE(i+1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;    
    
    fprintf('Finished iteration %d of %d in %4.1f s (MSE=%2.2e, cost=%2.2e)\n',i,iter,stats.time.iter(i+1),stats.MSE(i+1),stats.cost(i));
    if(stats.MSE(i+1) < settings.exit_tol)
    fprintf('**Reached exit tolerance: MSE < %2.2e\n',settings.exit_tol);
        break; 
    end
end
fprintf('Done! Total runtime: %6.1f s\n',sum(stats.time.iter));
fprintf('Final cost: %2.3e\n',stats.cost(end));
fprintf('Final MSE: %2.3e\n',stats.MSE(end));
fprintf('\n');
end



