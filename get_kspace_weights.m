function dz = get_kspace_weights(weighting,k,res,settings)
%predefined Fourier domain weighting schemes

if(strcmp(weighting,'grad')) %gradient weighting
    dz = zeros([res,2]);
    dz(:,:,1) = reshape(1j*2*pi*k(1,:),res)/max(res);
    dz(:,:,2) = reshape(1j*2*pi*k(2,:),res)/max(res);
elseif(strcmp(weighting,'deconv')); %inverse gaussian weighting
    a = 1/(2*settings.sig^2);
    kx = reshape(k(1,:),res);
    ky = reshape(k(2,:),res);
    gaussian = 2*sqrt(a/pi)*exp(-(pi^2)*(kx.^2 + ky.^2)/a);
    dz = 1./gaussian;     
elseif(strcmp(weighting,'none')) %no weighting
	dz = ones([res,1]);
else
    error('error: specify settings.weighting as grad, deconv, or none');
end

end