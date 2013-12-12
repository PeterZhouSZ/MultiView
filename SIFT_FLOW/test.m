clear all;
left = double(imread('left.bmp'));
frontal = double(imread('frontal.bmp'));
right = double(imread('right.bmp'));
[m,n] = size(left);
Z = zeros(1,m*n);
count = 1;
x = 1:m;
y = 1:n;
[X Y] = meshgrid(x,y);
% for i = 1:m
%     for j = 1:n
% %         Z(1,count) = i;
% %         Z(2,count) = j;
%         Z(i,j) = left(i,j);
%     end
% end
surface(X,Y,left')
figure(2)
surface(X,Y,right')
% options.ReducedDim = 10;
% [eigvector,eigvalue] = PCA(Z,options);
% Y = Z*eigvector;