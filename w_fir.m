%=========================================================================%

% IMAGE FILTERING CODE   : 1.0 

% Developers             : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform basic image filtering (FIR Wiener)
% Contact information    : moswal@stanford.edu

% Input options


% USE imabsdiff

% To check fot the correct number of input arguments in your function
% USE error(nargchk(low,high,nargin))


%=========================================================================%

% Clear previous memory items
clear all; 
close all;
clc;

disp('Extracting image and initializing variables')
tic;

% Grab the BLURRED image from the respective directory
im = imread...
    ('/Users/mehuloswal/me354_final_project/image_cereals/image2.jpg');

% Convert the image into a gray scale
im = rgb2gray(im);

% Scale the image by intensities [0,1]
im = mat2gray(im);

% Convert the image to double
im = im2double(im);

% Image size
[m,n] = size(im);

% Window size (M)                       [Later make this as an input]
M = 4;
W = -M:M;

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
% sig_n = 0.01;
sig_n = 0.0*sig;

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

% Covariance function (r): size = [2p+1,2p+1]
p = 3;                               % Need to check/play on this parameter
mpr = -p:p; npr = -p:p;

disp('Done')
toc;

%% Evaluating the RIGHT HAND SIDE OF Eqn. 8.64
disp('Evaluating the RIGHT HAND SIDE OF Eqn. 8.64')
tic;
% Creating h (PSF) AND r0 (Covariance)
r0 = zeros(length(mpr),length(npr));
h  = zeros(length(mpr),length(npr));

for i = -p:p
    for j = -p:p        
        r0(i+p+1,j+p+1) = 0.95^sqrt(i^2 + j^2);

        % PSF (h) : size = [2p+1,2p+1] 
        % Gaussian                      [Later make this as an input]
        alpha = .01; % This should be made smaller
        h(i+p+1,j+p+1) = exp(-alpha*(i^2 + j^2));
    end
end

% To use MATLAB's inbuilt function to create the PSF
% h = fspecial('disk',p);

% r_uv = Convolution between PSF and r0
r_uv = zeros(2*M+1,2*M+1);              % Initializing


for k = 1:length(W)
    for l = 1:length(W)
        for i = 0:length(mpr)-1
            for j = 0:length(npr)-1
                % Creating temp variables a and b
                a = k-i ; b = l-j;
                if a<=0 || a>length(mpr) || b<=0 || b>length(npr)
                    r_uv(k,l) = r_uv(k,l);
                else
                    r_uv(k,l) = r_uv(k,l) + h(a,b).*r0(i+1,j+1);
                end
            end         
        end        
    end
end

% Reshaping matrix r_uv into a vector
r_uv = reshape(r_uv,(2*M+1)^2,1);

disp('Done')
toc;

%% Evaluating the LEFT HAND SIDE OF Eqn. 8.64

disp('Evaluating the LEFT HAND SIDE OF Eqn. 8.64')
tic;
R = zeros((2*M+1)^2,(2*M+1)^2);
I = eye(size(R,1),size(R,2));

% Before we create R, we need to create matrix A (Autocorrelation of A)

A = zeros(2*M+1,2*M+1);
% Autocorrelation
for k = 1:length(W)
    for l = 1:length(W)
        for i = 0:length(mpr)-1
            for j = 0:length(npr)-1
                % Creating temp variables a and b
                a = k+i ; b = l+j;
                if a<=0 || a>length(mpr) || b<=0 || b>length(npr)
                    A(k,l) = A(k,l);
                else
                    A(k,l) = A(k,l) + h(a,b).*h(i+1,j+1);
                end
            end         
        end        
    end
end

% Now calculating the basic dimension of R

R_basic = zeros(length(W),length(W));
% Convolution
for k = 1:length(W)
    for l = 1:length(W)
        for i = 0:length(mpr)-1
            for j = 0:length(npr)-1
                % Creating temp variables a and b
                a = k-i ; b = l-j;
                if a<=0 || a>length(mpr) || b<=0 || b>length(npr)
                    R_basic(k,l) = R_basic(k,l);
                else
                    R_basic(k,l) = R_basic(k,l) + r0(a,b).*A(i+1,j+1);
                end
            end         
        end        
    end
end

%===========================BEGIN-NEW=====================================%
delta = eye(size(R_basic,1),size(R_basic,2));
LHS_basic = (1/SNR).*delta + R_basic;
row = [LHS_basic(1,:),zeros(1,(2*M+1)^2-length(LHS_basic(1,:)))];
col = [LHS_basic(:,1);zeros((2*M+1)^2-length(LHS_basic(1,:)),1)];

LHS = toeplitz(col,row);
%===========================END-NEW=======================================%

%===========================BEGIN-OLD=======================================%
% % Now to create a block Toeplitz matrix R
% row = [R_basic(1,:),zeros(1,(2*M+1)^2-length(R_basic(1,:)))];
% col = [R_basic(:,1);zeros((2*M+1)^2-length(R_basic(1,:)),1)];
% 
% R = toeplitz(col,row);
% 
% % Lets call matrix LHS : size((2*M+1)^2,(2*M+1)^2))
% LHS = (1/SNR).*I + R;
%===========================END-OLD=======================================%

% To get g, we need to solve the system of equations LHS*g = r_uv
g = LHS\r_uv;

% reshaping g back to a matrix
g = reshape(g,2*M+1,2*M+1);

disp('Done')
toc;

%% DFFT of the image

disp('Convolving blurred image with FIR Wiener Filter')
tic;
im_sharp = zeros(size(im));
% Convolution
for k = 1:m
    for l = 1:n
        for i = 0:length(W)-1
            for j = 0:length(W)-1
                % Creating temp variables a and b
                a = k-i ; b = l-j;
                if a<=0 || a>m || b<=0 || b>n
                    im_sharp(k,l) = im_sharp(k,l);
                else
                    im_sharp(k,l) = im_sharp(k,l) + im(a,b).*g(i+1,j+1);
                end
            end         
        end        
    end
end

disp('Done')
toc;

%% Sharpedned image

figure(1)
imshow(im)
title('original image')

figure(2)
imshow(im_sharp./max(im_sharp(:)))
title('Sharpened Image')
