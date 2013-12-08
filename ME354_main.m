%=========================================================================%

% IMAGE PROCESSING SCRIPT

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : The main code used to obtain the PSF using 
%                          various filters and compares them to the base
%                          case.
% Contact information    : deman@stanford.edu
%                          moswal@stanford.edu

% Input options

% 1. PSF_function       = 'gauss' uses a gaussian approximation 
%                       = 'disk' uses a disk approximation
%
% 2. plot               = 'yes' or 'no'
% 3. draw_kernel        = 'yes' or 'no'


%=========================================================================%

clc
close all
clear all

% Input options
PSF_function        = 'gaussian';
sigma               = 1;  
plot                = 'yes';
draw_kernel         = 'yes';
ClearIm             = 'DSC_0517.jpg'; %'peppers.png';            % Name of the clear image
BlurredIm           = 'DSC_0518.jpg'; %'peppers_blur.png';      % Name of the blurred image
noise               = 10^(-4);                   % Noise value (experience)




% Load and read images from a local directory
% Gets the images and calculates pre-processes them
At = imread(ClearIm);
A = rgb2gray(At);
%A = resize(A,1000);
FA = fft2(A);

Bt = imread(BlurredIm);
B = rgb2gray(Bt);
%B = resize(B,1000);
FB = fft2(B);

FH = blurr_func(A,B);

% Determines the PSF that will be the initial guess for the filters
% [M,N] = size(A);
% H = fspecial(PSF_function,[M,N],sigma);
% FH = abs(fft2(H));

%FH = fft2(C);


FB2 = FH.*FA;
B2 = ifft2(FB2);
B2 = real(B2)./max(max(abs(B2)));


%surf(real(FH))
figure 
imshow(A)
figure
imshow(B)
figure
imshow(B2)
%grid off

% figure
% surf(H,'EdgeColor','none')
figure
surf(real(ifft2(FH)),'EdgeColor','none')
