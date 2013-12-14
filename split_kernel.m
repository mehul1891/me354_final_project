%=========================================================================%
% SPLIT KERNEL           : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Finds the optical kernel between two images, and
% in this program is used to calculate the kernels of the parts of a
% previously split image
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
function [h] = split_kernel(u,v)
U = fft2(u); V = fft2(v);
G = V./U;
H = 1./G;
h = (fftshift(ifft2(H)));