%=========================================================================%
% FILTER SCRIPT          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Takes a 2D projection of the optical kernel and
%                          plots it out for comparizon
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
%
%   u = clear image
%   v = blurred image
%
% OUTPUT OPTIONS
%
%   h = the 2D projection of the optical kernel
%=========================================================================%

function [h] = real_kernel_2D_projection(u,v,color)
% color = 'k';
% u = u_312;

% Converts the blurred and the clear images to the frequency domain 
U = fft2(u);
V = fft2(v);

% Calculates the optical kernel by the use of the real and blurred image
H = V./U;

% Brings the kernel back to the time domain
H_real = (fftshift(ifft2(H)));

% Obtains the dimmensions of H_real an chooses the middle point to create a
% cut project that line to a 2D plane
[k1,k2] = size(H_real);
h = H_real(:,floor(k2/2)+1);

% Checks to make sure that the max peak is projected in 2D space
if max(max(H_real))~= max(h)
    disp('The cut is not aligned with the max')
end

plot(h,color,'LineWidth',1.5)
