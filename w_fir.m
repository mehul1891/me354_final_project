%=========================================================================%

% IMAGE FILTERING CODE: 1.0 

% Developer              : Mehul Oswal
% Organization           : Stanford University
% Objective              : Perform basic image filtering (FIR Wiener)
% Contact information    : moswal@stanford.edu

% Input options


% USE imabsdiff

% To check fot the correct number of input arguments in your function
% USE error(nargchk(low,high,nargin))


%=========================================================================%

% Clear previous memory items
clear all; close all;clc

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

% Window size (W)                       [Later make this as an input]
M = 4;
W = -M:M;

% St.D. of noise                        [Later make this as an input]
sig_n = 0.01;

% St.D of the image
sig = std2(im);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

% Covariance function (r): size = [2p+1,2p+1]
p = 2;                               % Need to check/play on this parameter
mpr = -p:p; npr = -p:p;

%% CREATING h AND r0

% Initializing Covariance and Gaussian PSF
r0 = zeros(length(mpr),length(npr));
h  = zeros(length(mpr),length(npr));


for i = -p:p
    for j = -p:p        
        r0(i+p+1,j+p+1) = 0.95^sqrt(i^2 + j^2);

        % PSF (h) : size = [2p+1,2p+1] 
        % Gaussian                      [Later make this as an input]
        alpha = 0.001; % This should be made smaller
        h(i+p+1,j+p+1) = exp(-alpha*(i^2 + j^2));
    end
end

% FIR Weiner Filter kernel

% Right hand side of Eqn 8.64 (r_uv)
% r_uv = Convolution between PSF and r0
r_uv = zeros(2*M+1,2*M+1);              % Initializing


for k = 1:length(W)
    for l = 1:length(W)
        for i = 1:length(mpr)
            for j = 1:length(npr)
                % Creating temp variables
                a = k-i ; b = l-j;
                if a<=0 || a>length(mpr) || b<=0 || b>length(npr)
                    r_uv(k,l) = r_uv(k,l);
                else
                    r_uv(k,l) = r_uv(k,l) + h(a,b)*r0(i,j);
                end
            
            end
         
        end
        
    end
end

%% Initializing R
% R = zeros((2*M+1)^2,(2*M+1)^2);
% I = eye(size(R,1),size(R,2));

R = conv2(r0,xcorr2(h,h));
size(R)

% Lets call this matrix A
% A = (1/SNR).*I + R;

for k = 1:2*M+1
    for l = 1:2*M+1
        a(k,l) = xcorr2(h(k,l),h(k,l));
        R1(k,l) = conv2(r0(k,l),a(k,l));
    end
end



% To get g, we have to use a for loop







% DFFT of the image 