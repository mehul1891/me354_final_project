%=========================================================================%

% IMAGE PROCESSING CODE: 1.0 

% Developers             : Mehul Oswal & David Manosalvas
% Organization           : Stanford University
% Objective              : Perform basic image convolution
% Contact information    : deman@stanford.edu

% Input options

% 1. kernel_function     = 'ext' to use lpfilter.m
%                        = 'int' to use our developed Gauss convolution
%                        kernel
%
% 2. boundary_condition  = 'periodic' to use periodic boundary conditions
%                        = 'none' to not use any bundary condition
%
% 4. plot                = 'yes' or 'no'
% 5. draw_kernel         = 'yes' or 'no'


%=========================================================================%


clc
clear all
close all

% Input options
kernel_function     = 'ext';
plot                = 'yes';
draw_kernel         = 'no';
boundary_condition  = 'periodic';
sigma               = 5;


% Load and read images from a local directory
A = imread('peppers.png');
A_red(:,:) = A(:,:,1);
A_green(:,:) = A(:,:,2);
A_blue(:,:) = A(:,:,3);

% Analyzing the image matrix
[M,N,O] = size(A);

% Switching between Kernel calculation options
switch kernel_function
    case 'ext'
        GF = lpfilter('gauss',M,N,sigma);
    case 'int'
%         Preferable to work in the frequency domain directly as we do not
%         have the physical dimensions of the image (x,y).

%         Here, x,y,xo,yo are some measure of wave number and hence
%         frequency (\omega)
        F = 1;
        xo = M/2; % Mean in x
        yo = N/2; % Mean in y
        gauss2D = @(x,y) F.*exp(-((x-xo).^2 + (y-yo).^2.)...
            /(2.*sigma.^2));
        [X,Y] = meshgrid(1:M,1:N);
        G = gauss2D(X,Y);
        GF = G';
    otherwise
        error('Wrong "kernel_function" selected')
end

% Switching between boundary conditions
switch boundary_condition
    case 'periodic'
        [k1,k2] = size(GF);
        GFpad = zeros(size(A_red));
        GFpad([end+1-floor(k1/2):end,1:ceil(k1/2)], ...
        [end+1-floor(k2/2):end,1:ceil(k2/2)]) = GF;
        clear GF
        GF = GFpad;
    case 'none'
        GF = fft2(GF);
        
otherwise
        error('Wrong "bondary_condition" selected')
end


GF = repmat(GF,[1,1,O]);

% Fourier Transform of the image
FA = fft2(A);

% Convolution in the spatial domain is product in the frequency domain
H = GF.*FA;


% Working back in the spatial domain
B = real(ifft2(H));

% Normalizing the matrix to position it between the required bounds (0-1)
% Maximum values per color
B_redmax = max(max(B(:,:,1)));
B_bluemax = max(max(B(:,:,2)));
B_greenmax = max(max(B(:,:,3)));

% Normalizing each matrix
B(:,:,1) = B(:,:,1)./B_redmax;
B(:,:,2) = B(:,:,2)./B_bluemax;
B(:,:,3) = B(:,:,3)./B_greenmax;

% Plotting
switch plot
    case 'yes'
        figure()
        imshow(A)
        title('Original image')

        figure()
        imshow(B)
        title('Noisy image')
    
    case 'no'
        disp('No images plotted')
    otherwise
        error('Wrong "plot" option selected')
end

% if (strcmp(draw_kernel,'yes')==1) 
%         figure 
%         contourf(GF(:,:,1))
%         title('Kernel')
%         colorbar
%         
%         figure 
%         contourf(ifft2(GF(:,:,1)))
%         title('Real Kernel')
%         colorbar
%         
% else
% end
