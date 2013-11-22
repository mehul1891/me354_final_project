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

% 2. plot                = 'yes' or 'no'


%=========================================================================%


clc
clear
close all

% Input options
kernel_function = 'ext';
plot            = 'yes';

% Load and read images from a local directory
A = imread('peppers.png');

% Analyzing the image matrix
[M,N,O] = size(A);


% Switching between Kernel calculation options
switch kernel_function
    case 'ext'
        sig = 10;
        G = lpfilter('gauss',M,N,sig);
    case 'int'
        F = 1;
        xo = 0;
        yo = 0;
        sig_x = 100;
        sig_y = 100;
        gauss2D = @(x,y) F.*exp(-((x-xo).^2./(2.*sig_x.^2)+...
            (y-yo).^2./(2.*sig_y.^2)));
        [X,Y] = meshgrid(linspace(1,N,N),linspace(1,M,M));
        G = gauss2D(X,Y);
        G = fft2(G);
    otherwise
        error('Wrong "kernel_function" selected')
end

G = repmat(G,[1,1,O]);

% Fourier Transform of the image
FA = fft2(A);

% Convolution in the spatial domain is product in the frequency domain
H = FA.*G;

% Working back in the spatial domain
B = real(ifft2(H));


% Plotting
switch plot
    case 'yes'        
        figure()
        imagesc(A)
        title('Original image')

        figure()
        imshow(B)
        title('Noisy image')

        figure 
        contourf(X,Y,G(:,:,1))
        title('Kernel')
        colorbar
    case 'no'
        disp('No images plotted')
    otherwise
        error('Wrong "plot" option selected')
end
