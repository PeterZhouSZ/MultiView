% Searching the reference face examples. The normalized face dataset is
% also generated and saved in this function by setting flag=0 in the first run.
function [neighbor,neighbor2,m,n]=findneighbor_pca(img1,k,input,input2,flag) %INDEX,
% input='F:\TraingSet\CMU_Database\Testing\Interpolated\Frontal\'; % The path for interpolated LR training samples  %Interpoalted_Frontal
% input2='F:\TraingSet\CMU_Database\Testing\Right_22.5\'; % The path for corresponding HR training samples


% [mm,nn]=size(flist);
imspace=[];                % matrix used for storing all interpolated LR samples
imspace2=[];               % matrix used for storing all HR samples

h=fspecial('gaussian',[49 49],49/2/3);

h2=ones(29,29);
h2=h2/(size(h2,1)*size(h2,2));

ext = '.png';
% mask=0;
% For the first time, set flag=0 to prepare for the dataset. If the dataset has been prepared and saved, set flag=1 to load the dataset;
% mm=60;
if(flag==0)
for i = 1:27
%     if i==INDEX
%         continue
%     else
    str_input = [input sprintf('%02d',i) ext];%strcat(input,num2str(i),'.bmp');
    str_input2 = [input2 sprintf('%02d',i) ext];%strcat(input2,num2str(i),'.bmp');
    img=imread(str_input);
    img2=imread(str_input2);

    img=double(img);
    img2=double(img2);
    [m,n]=size(img);

    [lm,lv]=local_mean_var(img,25);   % use the same parameter value as the one for input image
    img=img-lm;
    img=img./lv;
    img2=img2-lm;
    img2=img2./lv;
    %}
    %%%%%%%%%%%%%%%
    %img=I.*img;
    %img2=I.*img2;
    %img2=img2-img;
    %img=img-imfilter(img,h);
    %%%%%%%%%%%%%%%
    imspace=[imspace,reshape(img,m*n,1)];
    imspace2=[imspace2,reshape(img2,m*n,1)];
%     end
end
% Save the prepared dataset
save 'My_DEMO_TIP_4' imspace imspace2  m n  %mask
end
if(flag==1)
    variable=load('My_DEMO_TIP_4.mat');
    imspace=variable.imspace;
    imspace2=variable.imspace2;
    % mask=variable.mask;
    m=variable.m;
    n=variable.n;
end
% tic;
if(size(img1,1)==0)
%     test=imread([test_path,'000001.bmp']);
else
    test=img1;
end

test=double(test);
test=reshape(test,m*n,1);
% mask=mask/mm;
% mask=reshape(mask,m*n*mag*mag,1);

error=zeros(1,size(imspace,2));
%{
G=load('C:\Documents and Settings\Administrator\My Documents\MATLAB\imgdist_G.mat');
G=G.G;
%}
%
psi=mean(imspace,2);
L=imspace-repmat(psi,[1,size(imspace,2)]);
test=test-psi;
%
% Constructing the eigenface space. Run this only for the first time. 
if(flag==0)
[V,D]=eig(L'*L);
diagonal=diag(D);
[tmp,ind]=sort(diagonal,'descend');

dim=10;               % The number of eigenvectors to be used

V=V(:,ind(1:dim));
D=diag(diagonal(ind(1:dim)));
diagonal=diag(D);
diagonal=diagonal.^-0.5;
D2=diag(diagonal);
%{
L=single(L);
V=single(V);
D2=single(D2);
%}
% clear imspace imspace2 D

E=L*V*D2;
imspacew=E'*L;

save PCA_subspace_TIP E imspacew;
end
%}
if(flag==1)
load PCA_subspace_TIP;
end

testw=E'*test;
tic;
%

for i=1:length(error)
%     error(i)=(test-imspace(:,i))'*(test-imspace(:,i));
    error(i)=norm((testw-imspacew(:,i)),2);
%     error(i)=norm((test-imspace(:,i)).*mask,1);
    % error(i)=imgdist(test,imspace(:,i),G,2,m,n);
end
[dst,idx]=sort(error,'ascend');
if (nargin<2)
    k=1;
end
neighbor=imspace(:,idx(2:k));         % If the leave one out scheme is used, set it as imspace(:,idx(2:k)). Otherwise set it as imspace(:,idx(1:k))
neighbor2=imspace2(:,idx(2:k));       % the same as above
toc;

