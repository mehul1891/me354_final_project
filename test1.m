clc
clear
close all

A = imread('image1.jpg');
[M,N,O] = size(A);
sig = 100;
G = lpfilter('gauss',M,N,sig);
% F = .1;
% xo = 0;
% yo = 0;
% sig_x = 1;
% sig_y = 1;
% gauss2D = @(x,y) F.*exp(-((x-xo).^2./(2.*sig_x.^2)+(y-yo).^2./(2.*sig_y.^2)));
% [X Y] = meshgrid(linspace(-5,5,400),linspace(-5,5,393));

% G_mid = gauss2D(X,Y);
G = repmat(G,[1,1,3]);

FA = fft2(A);
% FG = fft2(G);
H = G.*FA;

B = real(ifft2(H));

figure
imagesc(A)

figure
imshow(B)
