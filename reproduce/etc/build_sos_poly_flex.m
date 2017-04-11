function sos = build_sos_poly_flex(U,s,q,res,h)
filter_siz = size(h);
filter_siz2 = 2*filter_siz-[1,1];

ind_filter = find(h);
normfac = prod(res)/prod(filter_siz);

k2 = get_kspace_inds(filter_siz2);
ind_filter_pad = get_lowpass_inds(k2,filter_siz);

sos_small = zeros(filter_siz2);
cj = zeros(filter_siz2);
for j=1:length(s)
    filter = zeros(filter_siz);
    filter(ind_filter) = U(:,j);
    cj(ind_filter_pad) = ifftshift(filter);
    sos_small = sos_small + ((1/s(j))^q)*(abs(ifft2(cj)).^2);
end

kout = get_kspace_inds(res);
ind_filter_out = get_lowpass_inds(kout,filter_siz2);

soshat = zeros(res);
soshat_small = fft2(sos_small);
soshat(ind_filter_out) = soshat_small;
sos = ifft2(soshat)/normfac;
end

