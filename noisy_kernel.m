%=========================================================================%
% NOISY OPTICAL KERNEL EXTRACTION          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Splits the blurred and clear images into 9
% parts, obtaines the optical kernel for each using a de-convolution, and
% to eliminate the noise, it averages the optical kernels to obtain an
% approximation of the true kernel.

% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS

%=========================================================================%
function [h]=noisy_kernel(u,v)
%u = imread('DSC_0517.jpg');
%v = imread('DSC_0518.jpg');

[u1 u2 u3 u4 u5 u6 u7 u8 u9] = picture_split(u);
[v1 v2 v3 v4 v5 v6 v7 v8 v9] = picture_split(v);

h1 = split_kernel(u1,v1);
h2 = split_kernel(u2,v2);
h3 = split_kernel(u3,v3);
h4 = split_kernel(u4,v4);
h5 = split_kernel(u5,v5);
h6 = split_kernel(u6,v6);
h7 = split_kernel(u7,v7);
h8 = split_kernel(u8,v8);
h9 = split_kernel(u9,v9);

h = (h1+h2+h3+h4+h5+h6+h7+h8+h9)/9;









