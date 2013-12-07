%=========================================================================%

% IMAGE FILTERING CODE   : 1.0 

% Developers             : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform image filtering (FIR Wiener)
% Contact information    : moswal@stanford.edu

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
read_image = 'global'; % Reading the good image from the local storage
blur_type  = 'disk';% Can be 'gaussian','disk', 'motion',etc
plot_original = 'no'; % Can be 'yes' or 'no'
plot_blurred  = 'yes'; % Can be 'yes' or 'no'

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
        RADIUS = 3;
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
mean_noise = 0;
var_noise = 0.0001;
im = imnoise(blurred, 'gaussian', ...
                mean_noise, var_noise);
            
if strcmp(plot_blurred,'yes')
    figure, imshow(im)
    title('Simulate Blur and Noise')
else
end
toc;

%%
disp('Extracting image and initializing variables')
tic;

% Image size
[m,n] = size(im);

% Window (W)                            [Later make this as an input]
M = 4; W = -M:M;

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
sig_n = sqrt(var_noise);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

% Covariance function (r): size = [2p+1,2p+1]
p = RADIUS;                         % Need to check/play on this parameter
mpr = -p:p; npr = -p:p;             % Square Matrix 

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
%         alpha = .01; % This should be made smaller
%         h(i+p+1,j+p+1) = exp(-alpha*(i^2 + j^2));
    end
end

% To use MATLAB's inbuilt function to create the PSF
h = fspecial('disk',p);

% Normalizing PSF & Covariance to have area under surface = 1;
h = h./sum(h(:));r0 = r0./sum(r0(:));

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

%% Sharpened image

figure()
imshow(im)
title('original image')

figure()
imshow(im_sharp./max(im_sharp(:)))
title('Sharpened Image')

%%
%===================USING MATLAB'S IFR WIENER=============================%

estimated_nsr = noise_var / var(I(:));
wnr3 = deconvwnr(im, PSF, estimated_nsr);
figure, imshow(wnr3)
title('Restoration of Blurred, Noisy Image Using MATLAB inbuilt function');

%===================USING MATLAB'S IFR WIENER=============================%
