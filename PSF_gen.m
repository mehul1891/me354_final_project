%=========================================================================%
% PSF SCRIPT             : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Generate Point Spread Functions (PSF)
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% 
% PSF_type               : 'gaussian'
%                          'disk'
%                          'motion'
% 
% PSF_size               : size of the PSF
% 
% factor                 : 'global'
%                          'local'
% H                      : psf (as a matrix)
%=========================================================================%

function [H] = PSF_gen(PSF_type,PSF_size,factor)

global gauss_size_factor disk_size_factor motion_size_factor
global gaussian_sigma

% Locally set these factors. Can be made automated based on iteration
switch factor
    case 'global'
        g_s_f = gauss_size_factor;
        d_s_f = disk_size_factor;
        m_s_f = motion_size_factor;
        g_s   = gaussian_sigma;
    case 'local'
        g_s_f = 2;
        d_s_f = 1;
        m_s_f = 2;
        g_s   = 5;
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
            g_s_f*PSF_size,g_s);
    case 'disk'
        psf = fspecial('disk',...
            d_s_f*PSF_size);
    otherwise
        error('Blur type specified not yet set up in the code')
end
H = psf;
end