

% TRIAL THURSDAY DEC 5TH NIGHT
size(im)


im1=im-mean(im(:));
size(im1)


imshow(im1)
im1f=fftn(im1);
rvv=ifftn(conj(im1f).*im1f);
noise=var_noise.*eye(size(im,1),size(im,2));
max(noise(:))


min(noise(:))


lhs = rvv-noisel;

lhs = rvv-noise;
pad_PSF = size(im);
H = fftn(PSF,pad_PSF);
size(H)


lhsf=fftn(lhs);
ruv=lhsf./H;
G=ruv./fftn(rvv);
im_sharp=real(ifftn(G.*im1));
imshow(im_sharp)
