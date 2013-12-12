% This is the main routine for the From Local Pixel Structure to Global Image Super-resolution (LPS2GIS) algorithm

% Before running the algorithm, you should already have a number of HR and
% LR (interpolated version) facial image pairs. These facial images should
% be roughly aligned based on the eyes' coordinates. The HR and LR images
% are stored in two different folders, and the pair of HR and LR images
% should have the same filename.      

% The paths for the above mentioned two folders will be used in LPS2GIS.m
% (line 17,18) and findneighbor_pca.m (line 4,5). The alogrithm will make use
% of those facial images for training and testing in a leave-one-out manner.  

% Type "LPS2GIS" in matlab to run the algorithm. Remember to add the path
% for all files into the matlab path setting.  

clear all; close all;clc;
addpath('utilties');
test_path='X:/MultiView/YaleB_exteneted/right_45/interpolated/4by4/'; % The path  where the interpolated LR images are stored. The images are used as the reference for searching and warping
% test_path_LH='D:\HZ\neighbor embedding\Frontal_P1\';
test_path_hr='X:/MultiView/YaleB_exteneted/right_45/phase1/'; % The path where the corresponding HR ground-truth are stored. The ground-truth are used for evaluating the reconstruction errors and providing the subsampled observation data.  
imout_path='X:/MultiView/YaleB_exteneted/right_45/Results/Liu/4by4/'; % The path where the reconstructed images will be saved.

if ~exist(imout_path, 'dir')
    mkdir(imout_path);
end

flist=dir([test_path,'*.png']);  % load the list for all input images. 
[count,oo]=size(flist);         % count is the total number of input images.

error_total=[];                 % error_total is used to store the reconstruction errors
% erc=ones(count,300);

mag = 4;                     % magnification factor. This should be consistent with the training images provided. 


for kk = 2:count             % reconstruct the input one by one
    pp = 0;
for sigma2=35              % We can reconstruct the input with different sigma values by setting sigma2=[30 35 ...]
for KK = 9                   % We can reconstruct the input with different numbers of reference face examples by setting KK=[7 8 ...]. As leave-one-out methodology is adopted, the used real number of face examples is KK-1. Because the HR image of the input face is excluded. 
fname=flist(kk).name;

img1=imread([test_path,fname]);

ini=imread([test_path_hr,fname]);

hh=fspecial('gaussian',[11 11],2);
ini=imfilter(ini,hh);

% Normalization
img1=double(img1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask
% mask0=imread('C:\Documents and Settings\Administrator\My Documents\MATLAB\mask0.bmp');
% % mask0=double(mask0)/255;
% img1=img1.*mask0;
% ini=double(ini).*mask0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%  Illumination normalizaton using Xie's method 
%
[lm,lv]=local_mean_var(img1,25);
toc;
img1=img1-lm;
img1=img1./lv;
%}


% Find the reference face examples in the gallery set. The path for the
% gallery set is set in the function "findneighbor_pca". img2, img2hr are
% the found interpolated LR and HR reference face examples. m,n are the size of
% the HR image.
[img2,img2hr,m,n]=findneighbor_pca(img1,KK,test_path, test_path_hr,1);


% Denormalization
%
img1=img1.*lv;
img1=img1+lm;
img2=img2.*repmat(lv(:),[1,size(img2,2)]);
img2=img2+repmat(lm(:),[1,size(img2,2)]);
img2hr=img2hr.*repmat(lv(:),[1,size(img2,2)]);
img2hr=img2hr+repmat(lm(:),[1,size(img2,2)]);
%}

img2=uint8(img2);
img2hr=uint8(img2hr);

tmp=[];
tmphr=[];
ermap=ones(m,n,size(img2,2));
for i=1:size(img2,2)
    tmp(:,:,i)=reshape(img2(:,i),m,n);
    % tmphr(:,:,i)=reshape(img2hr(:,i),m*mag,n*mag);
    tmphr(:,:,i)=reshape(img2hr(:,i),m,n);
end


% Begin to warp the reference images;
iii=repmat(img1,[1 1 3]);
img1=iii;
img1=uint8(img1);
tic;
%
for p=1:size(tmp,3)
    %
    iii=tmp(:,:,p);    
    iii=repmat(iii,[1 1 3]);
    iii=uint8(iii);
    [u, v, ermap(:,:,p)] = optic_flow_me(img1, iii); % Generating the flow field [u, v]
        
    iii=tmphr(:,:,p);
    iii=repmat(iii,[1 1 3]);
    iii = uint8( mywarp_rgb( double(iii), u, v ) ) ; % Warp the reference face examples in accordance with the input face.
    iii=reshape(iii(:,:,1),m,n);
    tmphr(:,:,p)=iii;
    %  
end
tmphr=double(tmphr);
ermap=double(ermap);
toc;


% Begin to learn the local pixel structures;
dx=3; % define the half size of the square neighborhood window
dy=3;
sigma=sigma2;
w=zeros(m,n,2*dx+1,2*dy+1);
%
tic;
[B,C]=wcg(tmphr,sigma,2,1,ermap);  % learn the local pixel structures;

w(2:end-1,2:end-1,1,1)=B{1,1};
w(2:end-1,2:end-1,1,2)=B{1,2};
w(2:end-1,2:end-1,1,3)=B{1,3};
w(2:end-1,2:end-1,2,1)=B{2,1};
w(2:end-1,2:end-1,2,2)=B{2,2};
w(2:end-1,2:end-1,2,3)=B{2,3};
w(2:end-1,2:end-1,3,1)=B{3,1};
w(2:end-1,2:end-1,3,2)=B{3,2};
w(2:end-1,2:end-1,3,3)=B{3,3};
toc;
%}


% ITERATIVE RECONSTRUCTION
% img1 is the input image to be super-resolved
% The constraint: w and img1;
% 1st way: solve iteratively
%%
img1hr=double(img1(:,:,1)); % img1hr is the estimation of the target HR image, and is initialized as the interpolated LR input.
img1hr(1:mag:end,1:mag:end)=double(ini(1:mag:end,1:mag:end)); % Reset the observation pixels' values

iter=25;
sd2=0;
sd=0;
%
% w(1:mag:end,1:mag:end,:,:)=0;
% w(1:mag:end,1:mag:end,2,2)=1;
%
ww=cell(3,3);
ww{1,1}=w(:,:,1,1);
ww{1,1}(1:end-1,1:end-1)=ww{1,1}(2:end,2:end);
ww{1,2}=w(:,:,1,2);
ww{1,2}(1:end-1,1:end)=ww{1,2}(2:end,1:end);
ww{1,3}=w(:,:,1,3);
ww{1,3}(1:end-1,2:end)=ww{1,3}(2:end,1:end-1);
ww{2,1}=w(:,:,2,1);
ww{2,1}(1:end,1:end-1)=ww{2,1}(1:end,2:end);
ww{2,2}=w(:,:,2,2);
ww{2,3}=w(:,:,2,3);
ww{2,3}(1:end,2:end)=ww{2,3}(1:end,1:end-1);
ww{3,1}=w(:,:,3,1);
ww{3,1}(2:end,1:end-1)=ww{3,1}(1:end-1,2:end);
ww{3,2}=w(:,:,3,2);
ww{3,2}(2:end,1:end)=ww{3,2}(1:end-1,1:end);
ww{3,3}=w(:,:,3,3);
ww{3,3}(2:end,2:end)=ww{3,3}(1:end-1,1:end-1);
%
tic;
tol=.01;
for mm=1:iter
    sd2=sd;
    sd=0;
    ss=zeros(m-2,n-2);
    tt=ww{1,1}.*img1hr;
    ss=ss+tt(1:end-2,1:end-2);
    tt=ww{1,2}.*img1hr;
    ss=ss+tt(1:end-2,2:end-1);
    tt=ww{1,3}.*img1hr;
    ss=ss+tt(1:end-2,3:end);
    tt=ww{2,1}.*img1hr;
    ss=ss+tt(2:end-1,1:end-2);
    tt=ww{2,2}.*img1hr;
    ss=ss+tt(2:end-1,2:end-1);
    tt=ww{2,3}.*img1hr;
    ss=ss+tt(2:end-1,3:end);
    tt=ww{3,1}.*img1hr;
    ss=ss+tt(3:end,1:end-2);
    tt=ww{3,2}.*img1hr;
    ss=ss+tt(3:end,2:end-1);
    tt=ww{3,3}.*img1hr;
    ss=ss+tt(3:end,3:end);
    dd=ss-img1hr(2:end-1,2:end-1);
    
%     ddpad=padarray(dd,[1 1]);
%     dd11=ddpad(1:end-2,1:end-2);
%     dd12=ddpad(1:end-2,2:end-1);
%     dd13=ddpad(1:end-2,3:end);
%     dd21=ddpad(2:end-1,1:end-2);
%     dd22=ddpad(2:end-1,2:end-1);
%     dd23=ddpad(2:end-1,3:end);
%     dd31=ddpad(3:end,1:end-2);
%     dd32=ddpad(3:end,2:end-1);
%     dd33=ddpad(3:end,3:end);
%     
%     dd0=dd11.*w(1:end-2,1:end-2,3,3)+dd12.*w(1:end-2,2:end-1,3,2)+dd13.*w(1:end-2,3:end,3,1)+dd21.*w(2:end-1,1:end-2,2,3)+dd23.*w(2:end-1,3:end,2,1)+dd31.*w(3:end,1:end-2,1,3)+dd32.*w(3:end,2:end-1,1,2)+dd33.*w(3:end,3:end,1,1);%+dd22.*w(2:end-1,2:end-1,2,2);

    img1hr(2:end-1,2:end-1)=img1hr(2:end-1,2:end-1)+0.3*(dd); %-0.3*dd0;
    sd=sum(sum(abs(dd)));
    % erc(kk,mm)=sd;
        
% Range constraint
idi=img1hr>255;
img1hr(idi)=255;
idi=img1hr<0;
img1hr(idi)=0;
img1hr(1:mag:end,1:mag:end)=double(ini(1:mag:end,1:mag:end)); % Reset the observation points

if(abs(sd-sd2)<tol)
    break;
end
end
toc;
%%
%{
error(kk)=calculate_error(uint8(img1hr),'MSE',[test_path_hr,fname])
error2(kk)=calculate_error(uint8(img1hr),'ESMSE',[test_path_hr,fname])
error3(kk)=calculate_error(uint8(img1hr),'MSSIM',[test_path_hr,fname])
imwrite(uint8(img1hr),[imout_path,'sigma_',int2str(sigma2),'_',fname])
%}
%
pp=pp+1;
error(pp)=calculate_error(uint8(img1hr),'MSE',[test_path_hr,fname])
error2(pp)=calculate_error(uint8(img1hr),'ESMSE',[test_path_hr,fname])
error3(pp)=calculate_error(uint8(img1hr),'MSSIM',[test_path_hr,fname])
% imwrite(uint8(img1hr),[imout_path,'sigma_',int2str(sigma2),'_',fname])
imwrite(uint8(img1hr),[imout_path,fname])
%
% imshow(img1hr,[0 255]);
end
% save([imout_path,'sigma_5_50_',int2str(kk)],'error','error2','error3');
% error_total=[error_total;error,0,error2,0,error3];
% save([imout_path,'209-228sigma_',int2str(sigma2)],'error','error2','error3')
% clear error error2 error3
end
error_total=[error_total;error,0,error2,0,error3];
clear error error2 error3;
end

 save([imout_path,'error_total_TIP'],'error_total');