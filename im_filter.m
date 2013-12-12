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
% INPUT OPTIONS
%=========================================================================%


function [u,G] = im_filter(v,filter_type,PSF_type,PSF_size,var_n,factor)

global gauss_size_factor disk_size_factor motion_size_factor

% Image variance
var_s = var(v(:));

% Signal to Noise ratio (SNR)
SNR = var_s/var_n;

% Locally set these factors. Can be made automated based on iteration
switch factor
    case 'global'
        g_s_f = gauss_size_factor;
        d_s_f = disk_size_factor;
        m_s_f = motion_size_factor;
    case 'local'
        g_s_f = 2;
        d_s_f = 1;
        m_s_f = 2;
    otherwise
        error('Wrong "factor" type specified')
end
        

% Internal PSF kernel to deblur the image with
switch PSF_type
    case 'motion'
        psf = fspecial('motion', ...
            m_s_f*PSF_size,...
            m_s_f*PSF_size);
    case 'gaussian'
        psf = fspecial('gaussian',...
            g_s_f*PSF_size,...
            g_s_f*PSF_size);
    case 'disk'
        psf = fspecial('disk',...
            d_s_f*PSF_size);
    otherwise
        error('Blur type specified not yet set up in the code')
end

% Making the PSF of the same size as the image
pad = size(v);
H = psf2otf(psf,pad);

% DFT of the image
V = fft2(v);

switch filter_type
    case 'inverse'
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
        s = 0.5;
        G = (H_inv.^s).*(conj(H)./(abs(H).^2 + 1/SNR)).^(1-s);
        U = G.*V;
        u = ifft2(U);
    case 'least_squares'
    case 'ED+filt'
        
        % This is just a trial method of sharpening. Want to see what edge 
        % + general filter does to G. Does it move closer to reality.
        
        disp('Enter the type of regular filter to use for this in quotes')
        gen_filt = input('');
        switch gen_filt
            case 'inverse'
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
end