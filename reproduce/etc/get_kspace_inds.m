function k = get_kspace_inds( res )
%GET_KSPACE_INDS define k-space index sets based on fft indexing
%input:  res = fourier domain resoultion
%output:   k = 2xprod(res) linear array of kspace indicies
if mod(res(2),2)
    indx = ifftshift(-((res(2)-1)/2):((res(2)-1)/2));
else
    indx = [0:((res(2)/2)-1), -(res(2)/2):-1];
end
if mod(res(1),2) 
    indy = ifftshift(-((res(1)-1)/2):((res(1)-1)/2));
else
    indy = [0:((res(1)/2)-1), -(res(1)/2):-1];
end
[kx,ky] = meshgrid(indx,indy);
k(1,:) = kx(:);
k(2,:) = ky(:);
end

