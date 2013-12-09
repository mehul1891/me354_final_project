%=========================================================================%

% IMAGE FILTERING CODE   : 1.4 

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
blur_type  = 'gaussian';% Can be 'gaussian','disk', 'motion',etc
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
M = 15; W = -M:M;
kmin=-(M-1);kmax=M-1;
lmin=-(M-1);lmax=M-1;

% St.D of the image
sig = std2(im);

% St.D. of noise                        [Later make this as an input]
sig_n = sqrt(var_noise);

% Signal to Noise ratio (SNR)
SNR = sig.^2/sig_n.^2;

% Covariance function (r): size = [2p+1,2p+1]
p = 4;                         % Need to check/play on this parameter
mpr = -p:p; npr = -p:p;             % Square Matrix 

disp('Done')
toc;

%% Evaluating the RIGHT HAND SIDE OF Eqn. 8.64
disp('Evaluating the RIGHT HAND SIDE OF Eqn. 8.64')
tic;


% for i = -p:p
%     for j = -p:p        
%         r0(i+p+1,j+p+1) = 0.95^sqrt(i^2 + j^2);
% 
%         % PSF (h) : size = [2p+1,2p+1] 
%         % Gaussian                      [Later make this as an input]
%          % This should be made smaller
%         h(i+p+1,j+p+1) = exp(-alpha*(i^2 + j^2));
%     end
% end
alpha = .001;


% r_uv = Convolution between PSF and r0
r_uv = zeros(2*M+1,2*M+1);              % Initializing


for k = -M:M
    for l = -M:M
        for i = k-p:k+p
            for j = l-p:l+p
                kk = k+M+1; ll = l+M+1;
                % Creating temp variables a and b
                a = k-i ; b = l-j;
                if a<-p || a>p || b<-p || b>p
                    r_uv(kk,ll) = r_uv(kk,ll);
                else
                    r_uv(kk,ll) = r_uv(kk,ll) + exp(-alpha*(a^2 + b^2))...
                        *.6^sqrt(i^2 + j^2);
                end
            end         
        end        
    end
end

% r_uv
% 
% r_uv2 = zeros(2*M+1,2*M+1); 
% 
% for k = -M:M
%     for l = -M:M
%         for i = -100:100
%             for j = -100:100
%                 kk = k+M+1; ll = l+M+1;
%                 % Creating temp variables a and b
%                 a = k-i ; b = l-j;
%                 if a<-p || a>p || b<-p || b>p
%                     r_uv2(kk,ll) = r_uv2(kk,ll);
%                 else
%                     r_uv2(kk,ll) = r_uv2(kk,ll) + exp(-alpha*(a^2 + b^2))...
%                         *0.95^sqrt(i^2 + j^2);
%                 end
%             end         
%         end        
%     end
% end
% r_uv2

% Reshaping matrix r_uv into a vector
r_uv = reshape(r_uv,(2*M+1)^2,1);

disp('Done')
toc;
%% LEFT HAND SIDE
A = zeros(2*M+1,2*M+1);
for k = -M:M
    for l = -M:M
        for i = -k-p:-k+p
            for j = -l-p:-l+p
                kk = k+M+1; ll = l+M+1;
                % Creating temp variables a and b
                a = k+i ; b = l+j;
                if a<-p || a>p || b<-p || b>p
                    A(kk,ll) = A(kk,ll);
                else
                    A(kk,ll) = A(kk,ll) + exp(-alpha*(i^2 + j^2)).*...
                        exp(-alpha*(a^2 + b^2));
                end
            end         
        end        
    end
end

%% R_basic
% Now calculating the basic dimension of R

R_basic = zeros(length(W),length(W));
% Convolution
for k = -M:M
    for l = -M:M
        for i = -M:M
            for j = -M:M
                % Creating temp variables a and b
                a = k-i ; b = l-j;
                kk = k+M+1; ll = l+M+1;
                ii = i+M+1;jj=j+M+1;
                if ii<=0 || ii>2*M+1 || jj<=0 || jj>2*M+1
                    R_basic(kk,ll) = R_basic(kk,ll);
                else
                    R_basic(kk,ll) = R_basic(kk,ll) + .6^sqrt(a^2 + b^2).*A(ii,jj);
                end
            end         
        end        
    end
end
% R_basic
% 
% R_basic2 = zeros(length(W),length(W));
% % Convolution
% for k = -M:M
%     for l = -M:M
%         for i = -100:100
%             for j = -100:100
%                 % Creating temp variables a and b
%                 a = k-i ; b = l-j;
%                 kk = k+M+1; ll = l+M+1;
%                 ii = i+M+1;jj=j+M+1;
%                 if ii<=0 || ii>2*M+1 || jj<=0 || jj>2*M+1
%                     R_basic2(kk,ll) = R_basic2(kk,ll);
%                 else
%                     R_basic2(kk,ll) = R_basic2(kk,ll) + 0.95^sqrt(a^2 + b^2).*A(ii,jj);
%                 end
%             end         
%         end        
%     end
% end
% 
% R_basic-R_basic2

%===========================BEGIN-NEW=====================================%
delta = eye(size(R_basic,1),size(R_basic,2));
LHS_basic = (1/SNR).*delta + R_basic;
row = [LHS_basic(1,:),zeros(1,(2*M+1)^2-length(LHS_basic(1,:)))];
col = [LHS_basic(:,1);zeros((2*M+1)^2-length(LHS_basic(1,:)),1)];

LHS = toeplitz(col,row);
%===========================END-NEW=======================================%

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

