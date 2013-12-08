%=========================================================================%

% IMAGE FILTERING CODE   : 4.0 

% Developers             : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform image filtering (Non Linear Root filter)
% Contact information    : moswal@stanford.edu
% Acknowledgements       : Prof. Girod, Stanford University

% Input options

% USE imabsdiff

% To check fot the correct number of input arguments in your function
% USE error(nargchk(low,high,nargin))


%=========================================================================%

%% Fresh Start
clear all; 
close all;
clc;

%% Input Options
read_image    = 'global'; % Reading the good image from the local storage
blur_type     = 'gaussian';% Can be 'gaussian','disk', 'motion',etc
plot_original = 'no'; % Can be 'yes' or 'no'
plot_blurred  = 'yes'; % Can be 'yes' or 'no'
add_noise     = 'no'; % Can be 'yes' or 'no'
alpha         = 1; % Can be less than 1 or greater than 1

%% Extracting image and introducing noise into it
disp('Introducing noise in the image')
tic;

switch read_image
    case 'local'
        I = imread...
        ('/Users/mehuloswal/me354_final_project/image_cereals/image1.jpg');
    case 'global'
        I = imread('peppers.png');
    case 'other'
        % Ask for the path of the image on the local system
        alt_path = input('Please specify alternate path of the image in single quotes:');
        I = imread(alt_path);
    otherwise
        error('Wrong "read_image" choice specified')
end


% Prepare the image: 
I = mat2gray(rgb2gray(im2double(I)));
if strcmp(plot_original,'yes')
    figure, imshow(I);
    title('Original Image');
else
end


switch blur_type
    case 'motion'
%       LEN = input('LEN=');
%       THETA = input('THETA=');
        LEN = 51;
        THETA = 11;
        PSF = fspecial('motion', LEN, THETA);
        blurred = imfilter(I, PSF, 'conv', 'circular');
    case 'gaussian'
%       ROW = input('No. of Rows=');
%       COL = input('No. of Cols=');        
        RADIUS = 10;
        ROW = RADIUS;
        COL = RADIUS;
        PSF = fspecial('gaussian', ROW, COL);
        blurred = imfilter(I, PSF, 'conv', 'circular');
    case 'disk'
%       RADIUS = input('Radius of the disk');
        RADIUS = 2;
        PSF = fspecial('disk', RADIUS);
        blurred = imfilter(I, PSF, 'conv', 'circular');
    otherwise
        error('Blur type specified not yet set up in the code')
end

% Additive noise
switch add_noise
    case 'yes'
        mean_noise = 0;
        var_noise = 0.0001; % on a scale of 0-1
        im = imnoise(blurred, 'gaussian', ...
                        mean_noise, var_noise);
    case 'no'
        var_noise = 0.0001; % on a scale of 0-1
        im = blurred;
    otherwise
        error('Wrong "add_noise" input choice')
end
            
if strcmp(plot_blurred,'yes')
    figure, imshow(im)
    title('Simulate Blur and Noise')
else
end
toc;

%% DON'T NEED THIS HERE
% Image size
[m,n] = size(im);

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
sig_n = sqrt(var_noise);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

%% DFT of Image
im_dft = fftn(im);

%% FFT of u_hat
im_sharp_dft = (abs(im_dft).^alpha).*exp(1j.*angle(im_dft));

%% IDFT of u_hat
im_sharp = real(ifftn(im_sharp_dft));

%% Plot
figure()
imshow(im_sharp)
title('Sharpened Image')
