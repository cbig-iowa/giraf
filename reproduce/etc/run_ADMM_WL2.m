function [x,resvec] = run_ADMM_WL2(x0,mu,M,Mt,AtA,MtM,Atb,gam,lambda,iter,tol)
x = x0;
Mx = M(x);
L = zeros(size(Mx));
ndz = size(L,3);
resvec = zeros(1,iter);
for i = 1:iter
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
    
    resvec(i) = norm(residue(:))/norm(Y(:));
    if (iter > 10) && (resvec(i) < tol)
        return;
    end
end

end

