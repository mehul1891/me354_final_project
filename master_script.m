%=========================================================================%
% MASTER SCRIPT          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Image retrieval, sharpness metric, test metric
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% 
%                      
%
% 
%                        
%
% 
% 
%=========================================================================%

%% Clear previous cache
clear all; close all; clc;

global gauss_size_factor disk_size_factor motion_size_factor
global gaussian_sigma

%% User Specified inputs

% Additive noise parameters
add_noise           = 'yes';
mean_n              =  0   ;
var_n               =  1e-5 ;


%% Images that this code is trained and tested with

nimages = 4;

% Test image 1: peppers.png
im1 = imread('barrels.jpg');

% Test image 2: lina.jpg
% Copyright of the PlayBoy magazine. Free redistribution cautioned.
im2 = imread('lena.tiff'); 

% Test image 3: cameraman.tif
im3 = imread('cameraman.tif');

% Test image 4: ssphere.jpg
im4 = imread('ssphere.jpg');

%% Make the images gray scale (reduces computational workload and within
% scope of this project). Also intensities are normalized to between [0,1]

for i = 1:nimages
    if nimages > 4
        error(...
        'no. of images must be checked again. If not, comment this out')
    else
        switch i
            case 1                
                Im1=im1;%Im1 = mat2gray(rgb2gray(im2double(im1)));
                im_size1 = size(Im1);
                lower_dimension(1) = min(im_size1);
            case 2                
                Im2 = mat2gray(rgb2gray(im2double(im2)));
                im_size2 = size(Im2);
                lower_dimension(2) = min(im_size2);
            case 3                
                Im3 = mat2gray(im2double(im3));
                im_size3 = size(Im3);
                lower_dimension(3) = min(im_size3);
            case 4                
                Im4 = mat2gray(im2double(im4));
                im_size4 = size(Im4);
                lower_dimension(4) = min(im_size4);
            otherwise
                error(...
            'User needs to intervene here in converting images gray scale')
        end
    end
end
%% Size of the PSF as a percentage of the lower dimension
PSF_factor = 0.01; % 4% of the lower dimension
for i = 1:nimages
    PSF_size(i) = ceil(PSF_factor*lower_dimension(i));
end

%% Now that we have the blur sizes, lets get the blur PSF + noise set up

% 1     : Gaussian PSF
% 2     : Disk PSF
% 3     : Motion PSF

% Additional factors can be set to 1 for consistency amongst blur types
% These factors are introduced because PSF_size = 1% was not significatn
gauss_size_factor    = 2;
gaussian_sigma       = 5;
disk_size_factor     = 1;
motion_size_factor   = 5;

for i = 1:nimages
    if nimages > 4
        error(...
        'no. of images must be checked again. If not, comment this out')
    else
        switch i
            case 1
                PSF_11 = fspecial('gaussian', ...
                    gauss_size_factor*PSF_size(i),...
                    gaussian_sigma);
                blurred_11 = imfilter(Im1, PSF_11, 'conv', 'circular');
                
                PSF_12 = fspecial('disk', disk_size_factor*PSF_size(i));
                blurred_12 = imfilter(Im1, PSF_12, 'conv', 'circular');
                
                PSF_13 = fspecial('motion', ...
                    motion_size_factor*PSF_size(i),...
                    motion_size_factor*PSF_size(i));
                blurred_13 = imfilter(Im1, PSF_13, 'conv', 'circular');
                
                switch add_noise
                    case 'yes'
                        blurred_11 = imnoise(blurred_11, 'gaussian', ...
                        mean_n, var_n);
                        blurred_12 = imnoise(blurred_12, 'gaussian', ...
                        mean_n, var_n);
                        blurred_13 = imnoise(blurred_13, 'gaussian', ...
                        mean_n, var_n);
                    case 'no'
                        disp('Train images only blurred, no noise added')
                    otherwise
                        error('Wrong "add_noise" input choice')
                end
            case 2
                PSF_21 = fspecial('gaussian', ...
                    gauss_size_factor*PSF_size(i), ...
                    gaussian_sigma);
                blurred_21 = imfilter(Im2, PSF_21, 'conv', 'circular');
                
                PSF_22 = fspecial('disk', disk_size_factor*PSF_size(i));
                blurred_22 = imfilter(Im2, PSF_22, 'conv', 'circular');
                
                PSF_23 = fspecial('motion', ...
                    motion_size_factor*PSF_size(i), ...
                    motion_size_factor*PSF_size(i));
                blurred_23 = imfilter(Im2, PSF_23, 'conv', 'circular');
                
                switch add_noise
                    case 'yes'
                        blurred_21 = imnoise(blurred_21, 'gaussian', ...
                        mean_n, var_n);
                        blurred_22 = imnoise(blurred_22, 'gaussian', ...
                        mean_n, var_n);
                        blurred_23 = imnoise(blurred_23, 'gaussian', ...
                        mean_n, var_n);
                    case 'no'
                        disp('Train images only blurred, no noise added')
                    otherwise
                        error('Wrong "add_noise" input choice')
                end
                
            case 3
                PSF_31 = fspecial('gaussian', ...
                    gauss_size_factor*PSF_size(i), ...
                    gaussian_sigma);
                blurred_31 = imfilter(Im3, PSF_31, 'conv', 'circular');
                
                PSF_32 = fspecial('disk', disk_size_factor*PSF_size(i));
                blurred_32 = imfilter(Im3, PSF_32, 'conv', 'circular');
                
                PSF_33 = fspecial('motion', ...
                    motion_size_factor*PSF_size(i), ...
                    motion_size_factor*PSF_size(i));
                blurred_33 = imfilter(Im3, PSF_33, 'conv', 'circular'); 
                
                switch add_noise
                    case 'yes'
                        blurred_31 = imnoise(blurred_31, 'gaussian', ...
                        mean_n, var_n);
                        blurred_32 = imnoise(blurred_32, 'gaussian', ...
                        mean_n, var_n);
                        blurred_33 = imnoise(blurred_33, 'gaussian', ...
                        mean_n, var_n);
                    case 'no'
                        disp('Train images only blurred, no noise added')
                    otherwise
                        error('Wrong "add_noise" input choice')
                end
                
            case 4
                PSF_41 = fspecial('gaussian', ...
                    gauss_size_factor*PSF_size(i), ...
                    gaussian_sigma);
                blurred_41 = imfilter(Im4, PSF_41, 'conv', 'circular');
                
                PSF_42 = fspecial('disk', disk_size_factor*PSF_size(i));
                blurred_42 = imfilter(Im4, PSF_42, 'conv', 'circular');
                
                PSF_43 = fspecial('motion', ...
                    motion_size_factor*PSF_size(i), ...
                    motion_size_factor*PSF_size(i));
                blurred_43 = imfilter(Im4, PSF_43, 'conv', 'circular'); 
                
                switch add_noise
                    case 'yes'
                        blurred_41 = imnoise(blurred_41, 'gaussian', ...
                        mean_n, var_n);
                        blurred_42 = imnoise(blurred_42, 'gaussian', ...
                        mean_n, var_n);
                        blurred_43 = imnoise(blurred_43, 'gaussian', ...
                        mean_n, var_n);
                    case 'no'
                        disp('Train images only blurred, no noise added')
                    otherwise
                        error('Wrong "add_noise" input choice')
                end
                
            otherwise
                error(...
            'User needs to intervene here in setting up image PSFs')
        end
    end
end


%% Plotting
% Index read guide
% u_221: 2nd image, disk blur, filter index

% Filter index
    % 1: inverse
    % 2: pseudo_inverse
    % 3: wiener
    % 4: geo_mean
    % 5: least_squares
    % 6: ED+filt

% % For Peppers 
%     % Gaussian
% v               = blurred_22;
% PSF_type        = 'disk';
% PSF_dim         = PSF_size(2);
% factor          = 'global';
% psf             = PSF_gen(PSF_type,PSF_dim,factor);
% 
% filter_type     = 'inverse';
% [u_111,G_111]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'pseudo_inverse';
% [u_112,G_112]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'wiener';
% [u_113,G_113]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'geo_mean';
% [u_114,G_114]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'least_squares';
% [u_115,G_115]   = im_filter(v,filter_type,psf,var_n);

% figure, imshow(real(u_112))
%title(['Peppers blurred by ', PSF_type ' image recovered using ',filter_type])

% Plots the 2D projection of the kernel

% figure
% [h_real] = real_kernel_2D_projection(real(Im2),v,'-b');
% hold on
% [h_111] = kernel_filter_2D_projection(G_111,'-r');
% [h_112] = kernel_filter_2D_projection(G_112,'-g');
% [h_113] = kernel_filter_2D_projection(G_113,'-y');
% [h_114] = kernel_filter_2D_projection(G_114,'-c');
% [h_115] = kernel_filter_2D_projection(G_115,':k');
% 
% error_11(1) = norm(h_real-h_111,2);
% error_11(2) = norm(h_real-h_112,2);
% error_11(3) = norm(h_real-h_113,2);
% error_11(4) = norm(h_real-h_114,2);
% error_11(5) = norm(h_real-h_115,2);
% 
% title(['Optical kernel for Peppers'])
% xlabel('pixels')
% ylabel('Magnitude')
% legend('Real', ['inverse filter ' num2str(error_11(1))],...
%     ['pseudo-inverse filter ' num2str(error_11(2))],...
%     ['wiener filter ' num2str(error_11(3))],...
%     ['geo-mean filter ', num2str(error_11(4))],...
%     ['least-squares filter ' num2str(error_11(5))]);


% Relevant plots
% For Lena 
    % disk
% v               = blurred_22;
% PSF_type        = 'disk';
% PSF_dim         = PSF_size(2);
% factor          = 'global';
% psf             = PSF_gen(PSF_type,PSF_dim,factor);
% 
% filter_type     = 'inverse';
% [u_111,G_111]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'pseudo_inverse';
% [u_112,G_112]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'wiener';
% [u_113,G_113]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'geo_mean';
% [u_114,G_114]   = im_filter(v,filter_type,psf,var_n);
% filter_type     = 'least_squares';
% [u_115,G_115]   = im_filter(v,filter_type,psf,var_n);
% 
% figure
% [h_real] = real_kernel_2D_projection(real(Im2),v,':b');
% hold on
% % [h_111] = kernel_filter_2D_projection(G_111,'-r');
% [h_112] = kernel_filter_2D_projection(G_112,'-k');
% [h_113] = kernel_filter_2D_projection(G_113,'-r');
% [h_114] = kernel_filter_2D_projection(G_114,'-c');
% % [h_115] = kernel_filter_2D_projection(G_115,':k');
% 
% % error_11(1) = norm(h_real-h_111,2);
% error_11(2) = norm(h_real-h_112,2);
% error_11(3) = norm(h_real-h_113,2);
% error_11(4) = norm(h_real-h_114,2);
% % error_11(5) = norm(h_real-h_115,2);
% 
% title(['Optical kernel for Lena'])
% xlabel('pixels')
% ylabel('Magnitude')
% legend('Real', ...
%     ['pseudo-inverse filter ' num2str(error_11(2))],...
%     ['wiener filter ' num2str(error_11(3))],...
%     ['geo-mean filter ', num2str(error_11(4))]);%,...
% %     ['least-squares filter ' num2str(error_11(5))]);

%===== Filter comparison for a balistic sphere with motion blurr ========%
% PSF charactristics for the initial guess
v               = blurred_43;
PSF_type        = 'motion';
PSF_dim         = PSF_size(4);
factor          = 'global';
psf             = PSF_gen(PSF_type,PSF_dim,factor);

% Filtering
filter_type     = 'inverse';
[u_111,G_111]   = im_filter(v,filter_type,psf,var_n);
filter_type     = 'pseudo_inverse';
[u_112,G_112]   = im_filter(v,filter_type,psf,var_n);
filter_type     = 'wiener';
[u_113,G_113]   = im_filter(v,filter_type,psf,var_n);
filter_type     = 'geo_mean';
[u_114,G_114]   = im_filter(v,filter_type,psf,var_n);
filter_type     = 'least_squares';
[u_115,G_115]   = im_filter(v,filter_type,psf,var_n);
[mssim ssim_map] = ssim_index(Im4, Im4);
disp(mssim)

% Plotting the 2D projection of the optical kernel
figure
[h_real] = real_kernel_2D_projection(real(Im4),v,':b');
hold on
% [h_111] = kernel_filter_2D_projection(G_111,'-r');
[h_112] = kernel_filter_2D_projection(G_112,'-k');
[h_113] = kernel_filter_2D_projection(G_113,'-r');
[h_114] = kernel_filter_2D_projection(G_114,'-c');
% [h_115] = kernel_filter_2D_projection(G_115,':k');
axis([170 180 -.005 .055])

% Norm2 used to calculate the difference betweent the kernels
% error_11(1) = norm(h_real-h_111,2);
error_11(2) = norm(h_real-h_112,2);
error_11(3) = norm(h_real-h_113,2);
error_11(4) = norm(h_real-h_114,2);
% error_11(5) = norm(h_real-h_115,2);

% Title | axis | legend
title(['Optical kernel for supersonic flow around a sphere'])
xlabel('pixels')
ylabel('Magnitude')
legend('Real', ...
    'pseudo-inverse filter ' ,...
    'wiener filter ' ,...
    'geo-mean filter ');%,...
%     ['least-squares filter ' num2str(error_11(5))]);

% Automated metrics used to compare the effect of the filters

[grad_11(1), simind_11(1), jn_11(1)] = sharpness_metrics(u_111, Im4);
[grad_11(2), simind_11(2), jn_11(2)] = sharpness_metrics(u_112, Im4);
[grad_11(3), simind_11(3), jn_11(3)] = sharpness_metrics(u_113, Im4);
[grad_11(4), simind_11(4), jn_11(4)] = sharpness_metrics(u_114, Im4);
[grad_11(5), simind_11(5), jn_11(5)] = sharpness_metrics(u_115, Im4);
[val_g_11,ind_g_11] = min(grad_11);
[val_s_11,ind_s_11] = max(simind_11);
[val_j_11,ind_j_11] = min(jn_11);
