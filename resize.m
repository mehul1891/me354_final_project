function [An] = resize(A,n)
[M,N] = size(A);
An = A(floor(M/2)-floor(n/2):floor(M/2)+floor(n/2),...
    floor(N/2)-floor(n/2):floor(N/2)+floor(n/2));
