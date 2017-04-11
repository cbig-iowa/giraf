function R = build_gram_matrix(gradx,filter_siz,filter_siz2,ind_filter2)
% Build gram matrix R=S(x)^*S(x) using ffts and neighborhood operators
    gradx_ifft = ifft2(gradx);
    sos = fft2(sum(abs(gradx_ifft).^2,3));
    sos2 = fftshift(reshape(sos(ind_filter2),filter_siz2));
    R = im2col(sos2,filter_siz);
    R = rot90(R,-1);
end

