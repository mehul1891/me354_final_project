%=========================================================================%

% IMAGE FILTERING CODE   : 3.0 

% Developers             : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform image filtering (Geometric Mean Filter)
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
s             =  0.; % Takes value between 0 and 1

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

%%
% Image size
[m,n] = size(im);

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
sig_n = sqrt(var_noise);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;


%% Zero padding for the image

disp('ZERO PADDING STILL NEEDS TO BE DONE if at all you are going to do it')


%% FIR Wiener Filter based on zero mean 

% pad_PSF = size(im);
% H = fftn(PSF,pad_PSF);

% Another option is to use : Save it for later
H = psf2otf(PSF,size(im));

H_inv = zeros(size(H));
for i=1:size(H,1)
    for j=1:size(H,2)
        if H(i,j)>=eps
            H_inv(i,j)=1/H(i,j);
        else
            H_inv(i,j)=H_inv(i,j);
        end
    end
end
            

G_s = (H_inv.^s).*(conj(H)./(abs(H).^2 + 1/SNR)).^(1-s);

%% DFT of Image
im_dft = fftn(im);

%% FFT of u_hat
im_sharp_dft = G_s.*im_dft;

%% IDFT of u_hat
im_sharp = ifftn(im_sharp_dft);

%% Plot
figure()
imshow(im_sharp)
title('Sharpened Image')
