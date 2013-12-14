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
%=========================================================================%

function [h] = kernel_filter_2D_projection(G, color)

% Takes the inverse of the kernel and obtains the optical kernel in the
% freq. domain
H = 1./G;

% This portion is used to eliminate all the inf and replace them by very
% small numbers. This was done primarily to compensate for the correction
% in the pseudo_inverse filter
[xi,yj] = size(H);
for i = 1:xi
    for j = 1:yj
        if H(i,j) > 9999
            H(i,j) = 10^(-1);
        end
    end
end

% Transforms the optical kernel from the frequency to the time domain

H_real = abs(fftshift(ifft2(H)));

% Obtains the dimmensions of G_real an chooses the middle point to create a
% cut project that line to a 2D plane
[k1,k2] = size(H_real);
h = H_real(:,floor(k2/2));

% Checks to make sure that the max peak is projected in 2D space
if max(max(H_real))~= max(h)
    disp('The cut is not aligned with the max')
end

plot(h,color,'LineWidth',1.5)
axis([floor(k1/2)+1-50 , floor(k1/2)+1+50 , 0 , max(h)+.05*max(h)])







