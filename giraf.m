function [x,cost] = giraf(xinit,b,A,At,sampmask,param,settings)
%GIRAF run GIRAF algorithm
%   Input:
%       xinit: initialization of x (in Fourier domain)
%       b: vector of sampled entries
%       A: function handle for sampling operator
%       At: function handle for transpose sampling operator
%       sampmask: logical mask of sampling locations
%       param: struct of GIRAF-specific parameters 
%           (see specification of fields in code below)
%       settings: struct of gloabl problem settings
%           (see specification of fields in code below)
%       
%   Output:
%       x: reconstructed Fourier data
%       cost: values of cost function per iteration
%
%   GIRAF solves the optimization problem:
%   min_x ||Ax-b||_2^2 + lambda*||T(x)||_p^p
%   where A is a undersampling operator, and
%   T(x) is a Toeplitz-like matrix built from k-space data x
%   and ||.||_p is the Schatten p quasi-norm, with 0 <= p <= 1.
%
%   For more details see the paper:
%   A Fast Algorithm for Convolutional Structured Low-Rank Matrix Recovery
%   G. Ongie & M. Jacob, 2017. 
%   Pre-print available online: 
%   https://arxiv.org/abs/1609.07429
%
%  Greg Ongie 4/10/2017

p = settings.p; %Schatten p value
q = 1-(p/2);    
lambda = settings.lambda; %regularization parameter
filter_siz = settings.filter_siz; %rectangular filter dimensions
res = settings.res; %recon pixel resolution
filter_siz2 = 2*filter_siz - [1,1]; %squared filter dimensions

iter = param.iter; %number of iterations
eta = param.eta;   %epsilon decrease factor;
eps = param.eps0;  %initial epsilon
epsmin = param.epsmin; %minimum epsilon value
ADMM_iter = param.ADMM_iter; %number of inner ADMM iterations
ADMM_tol = param.ADMM_tol;   %exit tolerance for ADMM
if(isfield(param,'delta'))   %ADMM conditioning parameter
    delta = param.delta;
else
    delta = 100; %default
end

if(isfield(param,'eps0')) %intial epsilson
    eps0 = param.eps0;
else
    eps0 = 0; %auto-initialize option
end

if(isfield(param,'overres')) %oversampled reconstruction grid
    overres = param.overres;
else
    overres = res + 2*filter_siz; %default
end

%intialize variables, operators, and index sets
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

x = xinit;
x_pad = zeros(overres);
x_pad(ind_full) = x;
xold = x;

cost = [];
fprintf('Starting GIRAF (p=%d, lambda = %1.1e)\n',p,lambda);
tic;
for i=1:iter   
    %step 1: Compute sos annihilating polynomial
    gradx = M(x_pad);
    G = build_gram_matrix(gradx,filter_siz,filter_siz2,ind_filter2);
    [U,S] = eig(G);
    ev = abs(diag(S));
    if i==1 %initialze epsilon
        if eps0 > 0 
            eps = eps0;
        else
            eps = 0.001*max(ev); %auto-init eps
        end
    end
    mu = build_sos_poly(U,ev+eps,q,overres,filter_siz,filter_siz2,ind_filter,ind_filter2);
        
    %step 2: ADMM solution of least squares problem
    gam = max(mu(:))/delta;  %set ADMM parameter some fraction of sos-mask max value
    [x_pad,~] = run_ADMM_WL2(x_pad,mu,M,Mt,sampmask_pad,MtMmask,Atb_pad,gam,lambda,ADMM_iter,ADMM_tol);
    x_pad = reshape(x_pad,overres);
    x = reshape(x_pad(ind_full),res);
    
    %update epsilon
    eps = max(eps/eta,epsmin);
    
    %cost computations (of previous iterate)
    if p == 0
        shatten = 0.5*sum(log(ev+epsmin));
    else
        shatten = (1/p)*sum((ev+epsmin).^(p/2));
    end        
    diff = A(xold)-b;
    if(lambda == 0) %equality constrained
        cost(end+1) = shatten;
    else
        cost(end+1) = norm(diff(:)).^2 + lambda*shatten;          
    end 
        
    %check stopping condition
    if(norm(x(:)-xold(:))/norm(xold(:)) < settings.exit_tol)
        fprintf('**Reached exit tolerance: relerr < %2.2e\n',settings.exit_tol);
        break; 
    end
    xold = x;
    
    fprintf('Finished iteration %d of %d (cost=%2.2e)\n',i,iter,cost(i));
end
runtime = toc;
fprintf('Done! Total runtime: %6.1f s\n',runtime);
fprintf('Final cost: %2.3e\n',cost(end));
fprintf('\n');

%Function to build sum-of-squares annihilation weights
function mu = build_sos_poly(U,s,q,overres,filter_siz,filter_siz2,ind_filter,ind_filter2)
    normfac = prod(overres)/prod(filter_siz);
    mu_small = zeros(filter_siz2);
    for j=1:length(s)
        filter_pad = zeros(overres);
        filter_pad(ind_filter) = ifftshift(reshape(U(:,j),filter_siz));
        filter = reshape(filter_pad(ind_filter2),filter_siz2);
        mu_small = mu_small + ((1/s(j))^q)*(abs(ifft2(filter)).^2);
    end
    muhat_small = fft2(mu_small);
    muhat = zeros(overres);
    muhat(ind_filter2) = muhat_small;
    mu = ifft2(muhat)/normfac;
end

% Function to build gram matrix G=T(x)^*T(x) 
% using ffts and neighborhood operators
function G = build_gram_matrix(gradx,filter_siz,filter_siz2,ind_filter2)
    gradx_ifft = ifft2(gradx);
    sos = fft2(sum(abs(gradx_ifft).^2,3));
    sos2 = fftshift(reshape(sos(ind_filter2),filter_siz2));
    G = im2col(sos2,filter_siz);
    G = rot90(G,-1);
end

%ADMM weighted least squares solver
%solves ||Ax-b||_2^2 + lambda*||D^{1/2}Mx||_2^2
%where D = diag(mu)
function [x,resvec] = run_ADMM_WL2(x0,mu,M,Mt,AtA,MtM,Atb,gam,lambda,iter,tol)
    x = x0;
    Mx = M(x);
    L = zeros(size(Mx));
    ndz = size(L,3);
    resvec = zeros(1,iter);
    for ii = 1:iter
        % Y subprob 
        Z = gam*(Mx+L);
        muinv = repmat((mu + gam).^(-1),[1,1,ndz]);
        Y = fft2(muinv.*ifft2(Z));

        % x subprob
        if(lambda == 0) %equality constrained
            y = Mt(Y-L)./(MtM + (MtM==0));
            x = Atb + y.*(1-AtA);
        else %lagrange formulation
            x = (Atb + lambda*gam*Mt(Y-L))./(AtA + lambda*gam*MtM);        
        end    

        % L update
        Mx = M(x);
        residue = Mx-Y;
        L = L + residue;

        resvec(ii) = norm(residue(:))/norm(Y(:));
        if (iter > 10) && (resvec(ii) < tol)
            return;
        end
    end

end

end



