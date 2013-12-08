%=========================================================================%

% IMAGE FILTERING CODE   : 2.0 

% Developers             : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform image filtering (FIR Wiener)
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
read_image    = 'local'; % Reading the good image from the local storage
blur_type     = 'other';% Can be 'gaussian','disk', 'motion',etc
plot_original = 'no'; % Can be 'yes' or 'no'
plot_blurred  = 'yes'; % Can be 'yes' or 'no'
add_noise     = 'no'; % Can be 'yes' or 'no''
edge_det      = 'no'; % Can be 'yes' or 'no'
var_noise     = 0.000001; % on a scale of 0-1

%% Extracting image and introducing noise into it
disp('Introducing noise in the image')
tic;

switch read_image
    case 'local'
        I = imread...
            ('/Users/mehuloswal/Documents/Dropbox/ME354/presentation/BlurredImageUsed.jpg');
%         ('/Users/mehuloswal/me354_final_project/image_cereals/image1.jpg');
        
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

%%

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
        RADIUS = 5;
        PSF = fspecial('disk', RADIUS);
        blurred = imfilter(I, PSF, 'conv', 'circular');
    case 'other'
        blurred = I;
%         LEN = 51;
%         THETA = 11;
%         PSF = fspecial('motion', LEN, THETA);
%         RADIUS = 1;
%         PSF = fspecial('disk', RADIUS);
        RADIUS = 5;
        ROW = RADIUS;
        COL = RADIUS;
        PSF = fspecial('gaussian', ROW, COL);
    otherwise
        error('Blur type specified not yet set up in the code')
end

% Additive noise
switch add_noise
    case 'yes'
        mean_noise = 0;
        var_n = 0.0001; % on a scale of 0-1
        im = imnoise(blurred, 'gaussian', ...
                        mean_noise, var_n);
    case 'no'
        
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

%%
% Image size
[m,n] = size(im);

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
sig_n = sqrt(var_noise);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

%% BLIND DE-CONVOLUTION: when the PSF is not known

Svv=log10(abs(fftshift(fft2(im))).^2 );

% USE MATLAB INBUILT 


%% Zero padding for the image

disp('ZERO PADDING STILL NEEDS TO BE DONE if at all you are going to do it')


%% FIR Wiener Filter based on zero mean 
pad = size(im);
% H = fft2(PSF,pad(1),pad(2));

% Another option is to use : Save it for later
H = psf2otf(PSF,size(im));

G = conj(H)./(abs(H).^2 + 1/SNR);

%% DFT of Image
im_dft = fftn(im);

%% FFT of u_hat
im_sharp_dft = G.*im_dft;

%% IDFT of u_hat
im_sharp = ifftn(im_sharp_dft);

%% Plot
figure()
imshow(im_sharp)
title('Sharpened Image')

%% Edge Detection + u_hat'
switch edge_det
    case 'yes'
        threshold = var(im(:));
        image_edge = edge(im,threshold);
        im_sh_edge = im_sharp+image_edge;
        figure, imshow(im_sh_edge);
    otherwise
        image_edge = zeros(size(im));
end

%% Finally recovered image






