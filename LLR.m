%local linear regression
clear all;
addpath('utilties');
Lpath = 'X:/MultiView/YaleB_exteneted/left_45/interpolated/4by4/';
Hpath=  'X:/MultiView/YaleB_exteneted/frontal/phase1/';
img_out_path='X:\MultiView\YaleB_exteneted\right_45\Results\Our\phase1\4by4\';
% Ipath='I:\TraingSet\CMU_Database\Testing\Interpolated\Frontal_4\Patch\10_10\';% 

% load([Hpath,'Hseries']);   % Hseries and Lseries contain the HR and LR patches respectively. Each patch is represented as a column vector. All HR/LR training patches form a large matrix. 
% load([Lpath,'Lseries']);  
if ~exist(img_out_path, 'dir')
    mkdir(img_out_path);
end

% Some configration
% m,n should be the two dimension of imresized img
s1 = 10; %7 %5  %20
s2 = 10;

overlap1 = 5; %4 %6 %10
overlap2 = 5;
d1 = s1-overlap1;
d2 = s2-overlap2;


inputpatch=[];

outpatch=[];
% img_path='F:\TraingSet\CMU_Database\Testing\Initial_Estimate\'; % path for LR image
% img_path_hr='F:\TraingSet\CMU_Database\Testing\Frontal\';% path for HR image which is used to evaluate the reconstruction error

ext = '.png';
% flist=dir([img_path,'*.',extension]);
% [count,qq]=size(flist);
count = 27;%68

m = 100;
n = 100;

fprintf('Generating the data......\n');
[H, L] = generateData(Lpath, Hpath, s1, s2, overlap1, overlap2, 1);
numOfpatch = size(H,2)/count;

%% Initial phase
fprintf('Initial phase.........\n')
for kk = 1:count
    tic;
    L1 = L;
    H1 = H;
    
    strPatch = (kk-1) * numOfpatch + 1;
    endPatch = kk * numOfpatch;
    inputpatch = L1(:,strPatch : endPatch); 
    L1(:,strPatch : endPatch) = [];
    H1(:,strPatch : endPatch) = [];
% % kk=1;
%     tic
% % str_in = strcat(img_path,num2str(kk), extension);  
% str_out= strcat(img_out_path,num2str(kk), extension);  
% img=imread(str_in);

% img=double(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% str_patch=strcat(Lpath,num2str(kk),'.mat');
% load(str_patch);
% inputpatch=patch2;            % Take Care
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load database pay attension to notation
% % end
% % [L1 H1]=generate_database(kk);
% str_H=strcat(Hpath,'Hseries_',num2str(kk), '.mat'); 
% load(str_H);
% H1=H;
% 
% clear H;
% 
% str_L=strcat(Lpath,'Hseries_',num2str(kk),'.mat'); 
% load(str_L);
% L1=H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inputpatch=single(inputpatch);
    [leng_row leng_col]=size(inputpatch);
    [LENG_ROW LENG_COL]=size(L1);
    pathc_select_L=[];
    pathc_select_H=[];
    for i=1:leng_col
        for j=1:count - 1 %67 69
          patch_select_L(:,j)=L1(:,i+(j-1)*leng_col);
          patch_select_H(:,j)=H1(:,i+(j-1)*leng_col);
        end
%         feaData = [inputpatch(:,i) patch_select_L];
%         [Y W] = lle(feaData,12,3);
%         weights = W(:,2:end);
        %W=MyLLR2(inputpatch(:,i),patch_select_L);
         outpatch(:,i)=pca_w(inputpatch(:,i),patch_select_L,patch_select_H);
    %       L_COMP = orthcomp(patch_select_L);
    %       W=(L_COMP'*inputpatch(:,i))/255;
%          outpatch(:,i)=patch_select_H*weights;
    %     
    end
    % outpatch=real(outpatch);
    % W=(L1'*inputpatch);
    % outpatch=H*W;
    R=recover(outpatch,s1,s2,m,n,overlap1);
%     imshow(R,[])
    % imshow(uint8(R));
    imwrite(uint8(R),fullfile(img_out_path, [sprintf('%02d', kk) ext]));
    fprintf('Pocessed number N=%d\n',kk);

% R1=eigentrans_multi('D:\HZ\neighbor embedding\Testing_Left\51.bmp','D:\HZ\neighbor embedding\Frontal\','D:\HZ\neighbor embedding\Left_22.5\',56,48,50);
% I=imread('D:\HZ\neighbor embedding\Testing_Frontal\51.bmp');
% subplot(1,4,1),imshow(I,[0 255]),colormap(gray);
% subplot(1,4,2),imshow(R,[0 255]),colormap(gray);
% subplot(1,4,3),imshow(R1,[0 255]),colormap(gray);
% subplot(1,4,4),imshow(I3,[0 255]),colormap(gray);
toc
end
