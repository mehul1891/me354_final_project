


u = imread('DSC_0517.jpg');
v = imread('DSC_0518.jpg');


[u1 u2 u3 u4 u5 u6 u7 u8 u9] = picture_split(u);
[v1 v2 v3 v4 v5 v6 v7 v8 v9] = picture_split(v);

[h1]=noisy_kernel(u1,v1);
[h2]=noisy_kernel(u2,v2);
[h3]=noisy_kernel(u3,v3);
[h4]=noisy_kernel(u4,v4);
[h5]=noisy_kernel(u5,v5);
[h6]=noisy_kernel(u6,v6);
[h7]=noisy_kernel(u7,v7);
[h8]=noisy_kernel(u8,v8);
[h9]=noisy_kernel(u9,v9);

h = (h1+h2+h3+h4+h5+h6+h7+h8+h9)/9;

surf(h,'LineStyle','none')
%%
[k1,k2] = size(h);

[k1,k2] = size(h);
h_line = h(:,floor(k2/2)+1);
figure
plot(h_line)
axis([floor(k1/2)+1-50 , floor(k1/2)+1+50 , 0 , max(h)+.05*max(h)])

