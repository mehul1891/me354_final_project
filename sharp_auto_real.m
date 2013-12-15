%=========================================================================%
% AUTOMATED SHARP SCRIPT : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Automated image sharpening, sharpness metrics.
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% 
% Blurred image and initial guesses for the PSF shape and size.                     
%
% 
% Outputs optimum filter type, sharpened image, and sharpness metrics                       
%
% 
% 
%=========================================================================%
%%
clear all; close all; clc
I = mat2gray(rgb2gray(im2double(imread('DSCN6155.jpg'))));
figure, imshow(I)
%%
guess = 9; % disk radius

guess1 = 50; % gaussian size
guess2 = 7; % gaussian sigma

iter = 1; iter_limit=2;

while iter<=iter_limit
for i = 2:4;
    
    PSF1 = fspecial('disk',guess);
%     PSF2 = zeros(size(PSF1));

    PSF2 = fspecial('gaussian',guess1,guess2);
%     PSF1 = zeros(size(PSF2));

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
    noise_var = 1e-4;
    estimated_nsr = noise_var / var(I(:));
    % 
    psf = PSF;
    var_n = noise_var;
    v = I;
    
    switch i
        case 1
            filter_type     = 'inverse';
            [u_111,G_111]   = im_filter(v,filter_type,psf,var_n);
            [gradient(i), simind(i), just_noticeable(i)] = ...
                sharpness_metrics(u_111, 0);
        case 2
            filter_type     = 'pseudo_inverse';
            [u_112,G_112]   = im_filter(v,filter_type,psf,var_n);
            [gradient(i), simind(i), just_noticeable(i)] = ...
                sharpness_metrics(u_112, 0);
        case 3
            filter_type     = 'wiener';
            [u_113,G_113]   = im_filter(v,filter_type,psf,var_n);
            [gradient(i), simind(i), just_noticeable(i)] = ...
                sharpness_metrics(u_113, 0);
        case 4
            filter_type     = 'geo_mean';
            [u_114,G_114]   = im_filter(v,filter_type,psf,var_n);
            [gradient(i), simind(i), just_noticeable(i)] = ...
                sharpness_metrics(u_114, 0);
        case 5
            filter_type     = 'least_squares';
            [u_115,G_115]   = im_filter(v,filter_type,psf,var_n);
            [gradient(i), simind(i), just_noticeable(i)] = ...
                sharpness_metrics(u_115, 0);
        otherwise
    end
end

% grad = gradient(2:4);JN=just_noticeable(2:4);
% opt_fil_gra = find(grad == min(grad))+1;
% opt_fil_JN = find(JN==min(JN))+1;

opt_fil_gra = find(gradient == min(gradient));
opt_fil_JN = find(just_noticeable==min(just_noticeable));

% For convergence
grad_inverse(iter) = gradient(1);
grad_pseudoinverse(iter) = gradient(2);
grad_wiener(iter) = gradient(3);
grad_geo_mean(iter) = gradient(4);
grad_least_squares(iter) = gradient(5);

JN_inverse(iter) = just_noticeable(1);
JN_pseudoinverse(iter) = just_noticeable(2);
JN_wiener(iter) = just_noticeable(3);
JN_geo_mean(iter) = just_noticeable(4);
JN_least_squares(iter) = just_noticeable(5);


figure, plot(gradient,'*--'), title('Gradient based sharpness metric');
figure, plot(just_noticeable,'*--r'), title('JNBM based sharpness metric')

if opt_fil_gra==opt_fil_JN && opt_fil_gra>=3;
    next_step = 'proceed';
else
    next_step = 'break';
end


% next_step = 'proceed';

switch next_step
    case 'proceed'
        switch opt_fil_gra
            case 1
                filter_type     = 'inverse';
                disp(['Optimum filter for these images is ',filter_type])
%                 figure, imshow(real(u_111)), title(['Image restored using ',filter_type]);
            case 2
                filter_type     = 'pseudo_inverse';
                disp(['Optimum filter for these images is ',filter_type])
                figure, imshow(real(u_112)), title(['Image restored using ',filter_type]);
            case 3
                filter_type     = 'wiener';
                disp(['Optimum filter for these images is ',filter_type])
                figure, imshow(real(u_113)), title(['Image restored using ',filter_type]);
            case 4
                filter_type     = 'geo mean';
                disp(['Optimum filter for these images is ',filter_type])
                figure, imshow(real(u_114)), title(['Image restored using ',filter_type]);
            case 5
                filter_type     = 'least squares';
                disp(['Optimum filter for these images is ',filter_type])
                figure, imshow(real(u_115)), title(['Image restored using ',filter_type]);
            otherwise
        end
        
        iterate = input('Should the program continue if you are not satisfied?');
        
%         iterate = 'y';
        switch iterate
            case 'n'
                iter = iter_limit+1;
            case 'y'
                
                guess = guess + 1; %disk radius 
                
                guess1 = guess1 + 10; % gaussian size
                guess2 = guess2 + 1; % gaussian sigma
                
                iter = iter + 1;
            otherwise
                error('answer only in strings "y" or "n"')
        end
    case 'break'
        disp('JNBN and Gradient metric do not agreee. Change something')
        opt_fil_gra
        opt_fil_JN
        iter = iter_limit + 1;
    otherwise
end
end