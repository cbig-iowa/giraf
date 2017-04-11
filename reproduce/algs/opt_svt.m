function [x,stats] = opt_svt(x0,b,A,At,sampmask,param,settings)
%Solves min_x 0.5*||Ax-b||_2^2 + lambda*||T(x)||_*
%where A is the measurement operator, and
%T(x) is the block Toeplitz matrix built from k-space data x
%This code is meant for benchmarking only
fprintf('Starting SVT (lambda = %1.2e, beta = %1.2e)\n',param.lambda,param.beta);

lambda = param.lambda; %regularization parameter
iter = param.iter; %number of iterations
beta = param.beta; %ADMM step-size parameter

filter_siz = settings.filter_siz; 
res = settings.res;

k = get_kspace_inds(res);
dz = get_kspace_weights(settings.weighting,k,res,settings);

[T, Th] = define_toeplitz_operators(dz,res,filter_siz,0);
ThT = Th(T(ones(res)));
ThTdc = ThT; ThTdc(1) = 1;
AtA = At(A(ones(res)));

Atb = At(b);
x = Atb;
%Z = zeros(size(T(x)));
G = zeros(size(T(x))); %lagrange multiplier

stats.time.svd_step = zeros(1,iter+1);
stats.time.iter = zeros(1,iter+1);
stats.MSE = zeros(1,iter+1);
stats.cost = zeros(1,iter);

stats.time.svd_step(1) = 0;
stats.time.iter(1) = 0;
stats.MSE(1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;

for i=1:iter
    tic;
    
    %step 1: svd soft-shrinkage step
    Z = T(x) + G;
    tic;
    [U, S, V] = svd(Z,'econ');
    stats.time.svd_step(i+1) = toc;
    s = diag(S);
    s0 = max(s-1/beta,0);
    Z = U*diag(s0)*V';         
    
    %step 2: least-squares problem
    if(lambda == 0) %equality constrained
        y = Th(Z-G)./ThTdc;
        x = Atb + y.*(1-AtA);
    else %lagrange formulation
        x = (Atb + lambda*beta*(Th(Z-G)))./(AtA + lambda*beta*ThT);       
    end
    
    %Lagrange multiplier update
    G = G + (T(x) - Z);
    
    stats.time.iter(i+1) = toc;
    
    %cost/MSE computations
    stats.MSE(i+1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;
    
    if(settings.cost_computations)
        [~,S,~] = svd(T(x),'econ');
        nuclear = norm(diag(S),1);
        diff = A(x)-b;
        if(lambda == 0) %equality constrained
            stats.cost(i) = nuclear;
        else
            stats.cost(i) = 0.5*norm(diff(:)).^2 + lambda*nuclear;          
        end
    end
    
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


