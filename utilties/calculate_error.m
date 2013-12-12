% Compute the reconstruction error with different choices of image quality measure.
function error=calculate_error(imout,imin, string)
% imout is the input image to be evaluated
% string is the term specifying the evaluation method. It can be one of the
% follows: "MSE","MSSIM"(is actually the SSIM) and "ESMSE"
% arg3 speicfies the file path where the ground-truth is.

% if (nargin<3)   % If the path for the ground-truth is not provided, we manually set it here.
% path1='D:\Huyu\experiment5\H-res\T\';
% path2='D:\Huyu\experiment5\test\';
% flist=dir([path2,'*.bmp']);
% % flist=dir([path1,'*.bmp']);
% img1=imread([path1,flist(1).name]);
% elseif (nargin==3)
%     img1=imread(arg3);
% end
img1= imin;
img2=imout;

% img1 is original image
% img2 is input

s1=2; % The pixels close to the image borders are excluded in the computation of the reconstruction errors. 
s2=2; %

if(strcmp(string,'PSNR'))
delta=(double(img1(s1*3:end-s1*3,s2*3:end-s2*3))-double(img2(s1*3:end-s1*3,s2*3:end-s2*3))).^2;
[m,n]=size(delta);
error=sum(sum(delta))/m/n;
error = 10*log10(255*255/error);
end

if(strcmp(string,'MSSIM'))
%{
[m,n]=size(img1(s1*3:end-s1*3,s2*3:end-s2*3));
w=9;
d=2;
c1=0.0001; c2=0.0009;
k=0;
ssim=0;
for mi=1:d:m-w
    for ni=1:d:n-w
        patch1=img1(mi:mi+w-1,ni:ni+w-1);
        patch2=img2(mi:mi+w-1,ni:ni+w-1);
        ave1=sum(sum(patch1))/w/w;
        ave2=sum(sum(patch2))/w/w;
        dev1=sqrt(sum((patch1(:)-ave1).^2)/(w*w-1));
        dev2=sqrt(sum((patch2(:)-ave2).^2)/(w*w-1));
        dev12=sum(sum((patch1-ave1).*(patch2-ave2)))/(w*w-1);
        k=k+1;
        ssim=ssim+(2*ave1*ave2+c1)*(2*dev12+c2)/(ave1^2+ave2^2+c1)/(dev1^2+dev2^2+c2);
    end
end
mssim=ssim/k;
error=mssim;
%}
    % Parameter setting
    %
    K = [0.01 0.03];
    window = fspecial('gaussian', 11, 1.5); % previous 1.5
    L = 255;
    %}
    [error, ssim_map] = ssim_index(img1(s1*3:end-s1*3,s1*3:end-s1*3), img2(s1*3:end-s1*3,s1*3:end-s1*3), K, window, L);
end

if(strcmp(string,'ESMSE'))
% img1=img1(s1*3:end-s1*3,s2*3:end-s2*3);
% img2=img2(s1*3:end-s1*3,s2*3:end-s2*3);
% Edge stability measure
var=cell(5,1);
var{1}=1.19;
var{2}=1.44;
var{3}=1.68;
var{4}=2.0;
var{5}=2.38;
% imgs1=cell(5,1);
% imgs2=cell(5,1);
map1=cell(5,1);
map2=cell(5,1);
% sobel_h=[-1 0 1;-2 0 2;-1 0 1];
% sobel_v=[-1 -2 -1;0 0 0;1 2 1];
c1map=zeros(size(img1));
c2map=zeros(size(img2));
for i=1:5
    map1{i}=edge(img1,'canny',[],var{i});
    map2{i}=edge(img2,'canny',[],var{i});
    % c1map=c1map+double(map1{i});
    % c2map=c2map+double(map2{i});
end 
max1=0;
max2=0;
num1=0;
num2=0;
previous1=1;
previous2=1;
[m,n]=size(c1map);
for x=1:m
    for y=1:n
        max1=0;
        num1=0;
        previous1=1;
    %%%%%%%%%%%%%%%%%
    for i=1:5
    if(map1{i}(x,y)==previous1)
        num1=num1+1;
        previous1=map1{i}(x,y);
    else
        if(num1>max1)
            max1=num1;
        end
        num1=0;
        previous1=1;
    end
    end
    %%%%%%%%%%%%%%%%
    c1map(x,y)=max1;
    end
end
for x=1:m
    for y=1:n
        max2=0;
        num2=0;
        previous2=1;
    %%%%%%%%%%%%%%%%%
    for i=1:5
    if(map2{i}(x,y)==previous2)
        num2=num2+1;
        previous2=map2{i}(x,y);
    else
        if(num2>max2)
            max2=num2;
        end
        num2=0;
        previous2=1;
    end
    end
    %%%%%%%%%%%%%%%%
    c2map(x,y)=max2;
    end
end

c1map=c1map(s1*3:end-s1*3,s2*3:end-s2*3);
c2map=c2map(s1*3:end-s1*3,s2*3:end-s2*3);
s=sum(sum((c1map+c2map)>0));
esmse=sum(sum((c1map-c2map).^2))/s;
error=esmse;
end



