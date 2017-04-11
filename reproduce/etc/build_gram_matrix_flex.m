function G = build_gram_matrix_flex(gradx,h,res)
% Build gram matrix G=S(x)^*S(x) using ffts and neighborhood operators
filter_siz = size(h);
filter_siz2 = 2*filter_siz-[1,1];
ind_filter = find(h);
k = get_kspace_inds(res);
ind_filter2 = get_lowpass_inds(k,filter_siz2);

gradx_ifft = ifft2(gradx);
sos = fft2(sum(abs(gradx_ifft).^2,3));
sos2 = fftshift(reshape(sos(ind_filter2),filter_siz2));
R = im2col(sos2,filter_siz);
R = rot90(R,-1);
G = R(ind_filter,ind_filter);

end

% function R = build_gram_matrix(gradx,filter_siz,filter_siz2,ind_filter2)
% % Build gram matrix R=S(x)^*S(x) using ffts and neighborhood operators
%     gradx_ifft = ifft2(gradx);
%     sos = fft2(sum(conj(gradx_ifft).*gradx_ifft,3));
%     sos2 = fftshift(reshape(sos(ind_filter2),filter_siz2));
%     R = im2col(sos2,filter_siz);
%     R = rot90(R,-1);
% end
% 
