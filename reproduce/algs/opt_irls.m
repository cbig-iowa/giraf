function [x,stats] = opt_irls(x0,b,A,At,sampmask,param,settings)
%Solves min_x 0.5*||Ax-b||_2^2 + lambda*||T(x)||_p^p
%where A is the measurement operator, and
%T(x) is the block Toeplitz matrix built from k-space data x
%and ||.||_p is the Shatten p-norm, with 0 <= p <= 1.
%This code is meant for benchmarking only
fprintf('Starting IRLS (p=%d, lambda = %1.2e)\n',param.p,param.lambda);

%lambda = param.lambda; %regularization parameter
iter = param.iter; %number of iterations
eta = param.eta; %epsilon decrease factor;
eps = param.eps;
epsmin = param.epsmin;
p = param.p;
q = 1-(p/2);   
lsqr_iter = param.lsqr_iter; %max LSQR iterations
lsqr_tol = param.lsqr_tol; %LSQR tolerance
iter = param.iter; %irls iterations

filter_siz = settings.filter_siz; 
res = settings.res;

k = get_kspace_inds(res);
dz = get_kspace_weights(settings.weighting,k,res,settings);

[T, Th] = define_toeplitz_operators(dz,res,filter_siz,0);

Atb = At(b);
Tb = T(Atb);
x = Atb;

ind_samples_comp = find(sampmask==0);
z = zeros(size(ind_samples_comp));

stats.time.svd_step = zeros(1,iter+1);
stats.time.iter = zeros(1,iter+1);
stats.MSE = zeros(1,iter+1);
stats.cost = zeros(1,iter);

stats.time.svd_step(1) = 0;
stats.time.iter(1) = 0;
stats.MSE(1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;

for i=1:iter
    tic;
    
    %step 1: weight matrix update
    tic;
    [~, S, V] = svd(T(x),'econ');
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
    W = V*(diag((s+eps).^(-q)))*V';
    
    %step 2: weighted least squares problem
    R = cglsR(T,Th,W,ind_samples_comp,res,size(Tb));
    B = -vec(Tb*W);
    z = cgls(R,B,0,lsqr_tol,lsqr_iter,[],z);
    x = zeros(res);
    x(ind_samples_comp) = z;
    x = x + Atb;

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

%function handle for CGLS
function R = cglsR(T,Th,W,ind_samples_comp,res,dims)
    vec = @(x) x(:);
    function y = afun(x,transp_flag)
       if transp_flag==2      % y = A'*x
          z = vec(Th(reshape(x,dims)*W'));
          y = z(ind_samples_comp);
       elseif transp_flag==1 % y = A*x
          z = zeros(res);
          z(ind_samples_comp) = x;
          y = vec(T(z)*W);
       end
    end
    R = @(x,transp_flag) afun(x,transp_flag);
end

end


