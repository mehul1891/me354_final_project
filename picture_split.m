%=========================================================================%
% IMAGE SPLITTING          : ME354 FINAL PROJECT, AUT 2013
%=========================================================================%

%=========================================================================%
% REPOSITORY INFORMATION

% Developers             : David Manosalvas & Mehul Oswal
% Organization           : Stanford University
% Objective              : Splits an image into 9 smaller images
% Contact information    : deman@stanford.edu & moswal@stanford.edu
%=========================================================================%

%=========================================================================%
% INPUT OPTIONS
% u = image
%=========================================================================%
function [piece1 piece2 piece3 piece4 piece5 piece6 piece7 piece8 piece9] = picture_split(u)

[M,N] = size(u);
 x_part_size = floor(M/3);
 y_part_size = floor(N/3);
 
 piece1 = u([1:x_part_size],[1:y_part_size]);
 piece2 = u([1+x_part_size:2*x_part_size],[1:y_part_size]);
 piece3 = u([1+2*x_part_size:3*x_part_size],[1:y_part_size]);
 
 piece4 = u([1:x_part_size],[1+y_part_size:2*y_part_size]);
 piece5 = u([1+x_part_size:2*x_part_size],[1+y_part_size:2*y_part_size]);
 piece6 = u([1+2*x_part_size:3*x_part_size],[1+y_part_size:2*y_part_size]);

 piece7 = u([1:x_part_size],[1+2*y_part_size:3*y_part_size]);
 piece8 = u([1+x_part_size:2*x_part_size],[1+2*y_part_size:3*y_part_size]);
 piece9 = u([1+2*x_part_size:3*x_part_size],[1+2*y_part_size:3*y_part_size]);
end