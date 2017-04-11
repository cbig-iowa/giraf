function [T, Th] = define_toeplitz_operators(dz,res,filter_siz,conjugate_flag)
    
    ndz = size(dz,3);
    conv_siz = res-filter_siz+[1,1];
    nrows = prod(conv_siz);
    ncols = prod(filter_siz);
    if(conjugate_flag)
        zdim = 2*ndz;
    else
        zdim = ndz;
    end
    dims = [zdim*nrows,ncols];
    scalefac = 1/sqrt(prod(dims));
          
    function Tx = forwardT(x)
        Tx = zeros(dims);
        for i=1:ndz
            D = x.*dz(:,:,i); %weight fourier data
            I = rot90(fftshift(D),2); %center and flip weighted fourier data
            rowind = ((i-1)*nrows+1):(i*nrows);
            Tx(rowind,:) = im2col(I,filter_siz,'sliding').';                 
        end
        if(conjugate_flag)
        for i=1:ndz                                       
            D = x.*dz(:,:,i);
            I = conj(fftshift(D)); %no flip, apply conjugate
            rowind = ((i+ndz-1)*nrows+1):((i+ndz)*nrows);
            Tx(rowind,:) = im2col(I,filter_siz,'sliding').';
        end 
        end
        Tx = scalefac*Tx;
    end

    function ThY = adjointT(Y)
        ThY = zeros(res);
        for i=1:ndz
            rowind = ((i-1)*nrows+1):(i*nrows);
            I = Y(rowind,:).';
            x = col2imstep(real(I),res,filter_siz,[1 1]) + 1i*col2imstep(imag(I),res,filter_siz,[1 1]);
            x = ifftshift(rot90(x,-2)).*conj(dz(:,:,i));
            ThY = ThY + x;                        
        end
        if(conjugate_flag)
        for i=1:ndz
            rowind = ((i+ndz-1)*nrows+1):((i+ndz)*nrows);
            I = Y(rowind,:).';
            x = col2im(I,filter_siz,res,'sliding');
            x = ifftshift(conj(x)).*conj(dz(:,:,i));
            ThY = ThY + x;                        
        end
        end        
        ThY = scalefac*ThY;
    end
          
    T = @(x) forwardT(x);
    Th = @(Y) adjointT(Y);

end

 