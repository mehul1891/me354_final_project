% Trial File to play with filters.

A = im2double(...
    rgb2gray(...
    imread('/Users/mehuloswal/me354_final_project/new_pictures/DSCN6155.JPG')));
% I = A(973:2050,2509:3605);
I=A;
figure, imshow(I)

% Two PSF's added on top of each other to get the best sharpness
PSF1 = fspecial('disk',10);
% PSF2 = zeros(size(PSF1));
PSF2 = fspecial('gaussian',60,8);


% making the PSFs of the same size. Can also be considered as zero
% padding.

nrows = max(size(PSF1,1), size(PSF2,1));
ncols = max(size(PSF1,2), size(PSF2,2));

row_cen = ceil(nrows./2);
col_cen = ceil(ncols./2);

nchannels = size(PSF1,3);

extendedPSF1 = zeros(nrows,ncols);
extendedPSF2 = zeros(nrows,ncols);

% Temporary variables
a = ceil(size(PSF1,1)./2);
b = size(PSF1,1)-a;
aa = ceil(size(PSF1,2)./2);
bb = size(PSF1,2)-aa;

c = ceil(size(PSF2,1)./2);
d = size(PSF2,1)-c;
cc = ceil(size(PSF2,2)./2);
dd = size(PSF2,2)-cc;


extendedPSF1(row_cen-a+1:row_cen+b,col_cen-aa+1:col_cen+bb) = PSF1; 
extendedPSF2(row_cen-c+1:row_cen+d,col_cen-cc+1:col_cen+dd) = PSF2; 

PSF1=extendedPSF1;
PSF2=extendedPSF2;


PSF = PSF1+PSF2;
noise_var = 0.001;
estimated_nsr = noise_var / var(I(:));


wnr2 = deconvwnr(I, PSF, estimated_nsr);

figure, imshow(wnr2)