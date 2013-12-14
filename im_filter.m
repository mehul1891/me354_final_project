%=========================================================================%
% FILTER SCRIPT          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Filter the image using different kernels
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT & OUTPUT OPTIONS
 
% Input image (v)        : Blurry+noisy image to be sharpened
% 
% filter_type            : 'inverse'
%                          'wiener'
%                          'geom_mean'
%                          'least_squares'
%                          'ED+filt'
% 
% var_n                  : Some estimate of the variance of noise.Usually
%                          between 1e-4 to 1e-8
% 
% Output image (u)       : Sharpened image output 
% 
% filter kernel (G)      : Filter sharpen kernel in fourier domain
%=========================================================================%


function [u,G] = im_filter(v,filter_type,psf,var_n)

% Image variance
var_s = var(v(:));

% Signal to Noise ratio (SNR)
SNR = var_s/var_n;

% Making the PSF of the same size as the image
pad = size(v);
H = psf2otf(psf,pad);

% DFT of the image
V = fft2(v);

switch filter_type
    case 'inverse'
        G = 1./H;
        U = G.*V;
        u = ifft2(U);
        
    case 'pseudo_inverse'
        [xi,yj] = size(H);
        for i = 1:xi
            for j = 1:yj
                if H(i,j) <= 10^(-7)
                    H(i,j) = inf;
                end
            end
        end
        G = 1./H;
        U = G.*V;
        u = ifft2(U);
        
    case 'wiener'
        G = conj(H)./(abs(H).^2 + 1/SNR);
        U = G.*V;
        u = ifft2(U);
    
    case 'geo_mean'
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
        s = 0.52;
        G = (H_inv.^s).*(conj(H)./(abs(H).^2 + 1/SNR)).^(1-s);
        U = G.*V;
        u = ifft2(U);
    
    case 'least_squares'
        NP = var_n*prod(size(v));
        u = deconvreg(v,psf,NP);
        G = fft2(u)./V;
    
    case 'ED+filt'
        
        % This is just a trial method of sharpening. Want to see what edge 
        % + general filter does to G. Does it move closer to reality.
        
        disp('Enter the type of regular filter to use for this in quotes')
        gen_filt = input('');
        
        switch gen_filt
            case 'inverse'
                G_filt = 1./H;
                U_filt = G_filt.*V;
                u_filt = ifft2(U_filt);

            case 'pseudo_inverse'
                [xi,yj] = size(H);
                for i = 1:xi
                    for j = 1:yj
                        if H(i,j) <= 10^(-1)
                            H(i,j) = inf;
                        end
                    end
                end
                G_filt = 1./H;
                U_filt = G_filt.*V;
                u_filt = ifft2(U_filt);
            case 'wiener'
                G_filt = conj(H)./(abs(H).^2 + 1/SNR);
                U_filt = G_filt.*V;
                u_filt = ifft2(U_filt);
            case 'geo_mean'
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
                s = 0.5;
                G_filt = (H_inv.^s).*(conj(H)./(abs(H).^2 + 1/SNR)).^(1-s);
                U_filt = G_filt.*V;
                u_filt = ifft2(U_filt);
            case 'least_squares'
                NP = var_n*prod(size(v));
                u_filt = deconvreg(v,psf,NP);
            otherwise
            error('Check type of filter selected')
        end
        
        % threshold = var(im(:));
        image_edge = edge(v);
        u_edge = u_filt+image_edge;
        G = fft2(u_edge)./V;
        u = u_edge;
    
    otherwise
        error('Wrong "filter_type" specified')
end

% Normalizes the cleared image
u = u./max(max(u));
end