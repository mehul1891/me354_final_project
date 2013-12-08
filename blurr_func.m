%=========================================================================%

% FINDS THE FUNCTION THAT DEFINES THE OPTICS 

% Developers             : David Manosalvas
% Organization           : Stanford University
% Objective              : This function performs a decombolution to find 
%                          a representation of the optics that caused 
%                          blurr.
% Contact information    : deman@stanford.edu

% Input options
% A - The clear picture that will be used as a reference
% B - Blurry picture

% Output options
% C_real - The function that represents the optics and has been 
%          convolved with the image to generate the blurry picture

%=========================================================================%
function [FC] = blurr_func(A,B)

% Takes only the red part of the picture and its based on the assumption
% that the optics function will affect all the layers equaly
x(:,:) = A(:,:,1);
y(:,:) = B(:,:,1);

FA = fft2(x);
FB = fft2(y);

% Linear decombolution
FC = FB./FA;

%C = ifft2(FC); % The periodic BC based combolution function
%[k1,k2] = size(C);
%C_real = abs(fftshift(ifft2(FC)));

% Obtaines the DC centered signal
% C_real = C([end+1-floor(k1/2):end,1:ceil(k1/2)], ...
%         [end+1-floor(k2/2):end,1:ceil(k2/2)]);




