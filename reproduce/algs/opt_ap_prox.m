function [x,stats] = opt_ap_prox(x0,b,A,At,sampmask,param,settings)
% AP-PROX: Alternating projection with proximal smoothing
% Solves min_x 0.5*||Ax-b||_2^2 + lambda*phi(T(x))
% where phi(Y) := min_{Z: rank Z <= r} ||Y-Z||_2^2 
% (i.e. the l2-norm of the residual of best rank r approximation),
% A is a Fourier undersampling mask, and
% T(x) is the block Toeplitz matrix built from k-space data x
fprintf('Starting AP-PROX (r = %d, lambda = %1.2e)\n',param.r,param.lambda);

%addpath('etc/lansvd'); %include lansvd package

r = param.r; %rank cutoff
lambda = param.lambda; %regularization parameter
iter = param.iter; %number of iterations

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

stats.time.svd_step = zeros(1,iter+1);
stats.time.iter = zeros(1,iter+1);
stats.MSE = zeros(1,iter+1);
stats.cost = zeros(1,iter);

stats.time.svd_step(1) = 0;
stats.time.iter(1) = 0;
stats.MSE(1) = (norm(x(:)-x0(:))/norm(x0(:)))^2;
for i=1:iter
    tic;
    
    %step 1: projection onto set of rank r matricies
    Tx = T(x);
    tic; 
    switch param.svd_type
        case 'econ'
            [U, S, V] = svd(Tx,'econ');
        case 'lansvd'
            [U, S, V] = lansvd(Tx,r,'L');
        case 'rsvd'
            [U, S, V] = rsvd(Tx,r);            
    end
    stats.time.svd_step(i+1) = toc;
    s = diag(S);
    ranknorm = 0.5*norm(s((r+1):end))^2;    
    s((r+1):end) = 0;
    Z = U*diag(s)*V';
    
    if(settings.cost_computations)
        diff = A(x)-b;
        if(lambda == 0) %equality constrained
            stats.cost(i) = ranknorm;
        else
            stats.cost(i) = 0.5*norm(diff(:)).^2 + lambda*ranknorm;          
        end    
    end
    
    %step 2: least-squares problem
    if(lambda == 0) %equality constrained
        y = Th(Z)./ThTdc;
        x = Atb + y.*(1-AtA);
    else %lagrange formulation
        x = (Atb + lambda*(Th(Z)))./(AtA + lambda*ThT);       
    end    
    
    stats.time.iter(i+1) = toc;
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


