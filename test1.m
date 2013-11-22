clc
clear
close all

%A = imread('image1.jpg');
A = imread('peppers.png');
[M,N,O] = size(A);
sig = 10;
G = lpfilter('gauss',M,N,sig);
% F = 1;
% xo = 0;
% yo = 0;
% sig_x = 100;
% sig_y = 100;
% gauss2D = @(x,y) F.*exp(-((x-xo).^2./(2.*sig_x.^2)+(y-yo).^2./(2.*sig_y.^2)));
[X Y] = meshgrid(linspace(1,N,N),linspace(1,M,M));
% 
% G = gauss2D(X,Y);
G = repmat(G,[1,1,3]);

FA = fft2(A);
%G = fft2(G);

H = FA.*G;

B = real(ifft2(H));

figure
imagesc(A)

figure
imshow(B)

figure 
contourf(X,Y,G(:,:,1))
colorbar
