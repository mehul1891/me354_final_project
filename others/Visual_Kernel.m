% Script used to create Kernels for visualization only

clc
clear
close all

x = -3:.1:3;
y = 1*exp(-x.^2);
y1 = .8*(y +.2.*sin(2*pi*x).*y);
y2 = y +.1*rand(1,length(x));
y3 = 1*(y +.2.*sin(4*pi*x).*y); 

plot(x,y,'LineWidth',2)
hold on
plot(x,y1,'r','LineWidth',2)
plot(x,y2,'g','LineWidth',2)
plot(x,y3,'k','LineWidth',2)
axis ([-3 3 0 1.2])
legend('True Kernel','Filter 1', 'Filter 2', 'Filter 3')
title('Optical Kernel Comparison')

