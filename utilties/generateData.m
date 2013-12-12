% Designed for SR using neighbor embedding [1]	C. Hong, Y. Dit-Yan, and X. Yimin,¡°Super-resolution through neighbor embedding,¡± in Proc. IEEE Conf. Computer Vision and Pattern Recognition, 2004, pp. I-275-I-282." 
% function epatch()
% The path 'input' is for low-resolution, input2 is for high-res 
% The 'patch' represents the low-res while patch2 represents high-res
function [H, L] = generateData(input, input2, s1,s2,overlap1,overlap2,mag)

% input='E:\TraingSet\CMU_Database\Testing_Interpolated_Frontal\';
% input2='E:\TraingSet\CMU_Database\Testing_Frontal\';

%
% output='E:\TraingSet\CMU_Database\Testing_Interpolated_Frontal\Patch_NE\';   %LR
% output2='E:\TraingSet\CMU_Database\Testing_Frontal\Patch_NE\';   %HR
ext='png';
flist=dir([input,'*.',ext]);
[mm, nn]=size(flist);
% mag=1;
% 
% s1=7;
% s2=7;
ss1=s1*mag;
ss2=s2*mag;
% overlap1=6;
% overlap2=6;
d1=s1-overlap1;
d2=s2-overlap2;
dd1=d1*mag;
dd2=d2*mag;
k=0;

L = [];
H = [];
% h2= fspecial('gaussian',[15,15],18/2/3);
for i=1:mm
    patch1=[];
    patch2=[];
    img = imread([input sprintf('%02d.', i) ext]);
    img2=imread([input2 sprintf('%02d.', i) ext]);
    %img=imread([input num2str(i) ext]);
    %img2=imread([input2,num2str(i) ext]);
    img=double(img);
    img2=double(img2);
    %img2=imresize(img2,size(img),'nearest');
    [m,n]=size(img);
    % img3=imread([input3,flist(i).name(1:end-4),'.bmp']);
    % bw=edge(img3,'canny',[0.001, 0.0012]);
    % bw=ones(size(img3));
    % x=1:d1:m;
    % y=1:d2:n;
    % patch=zeros(s1*s2,length(y)*length(x));
    k=0;
    %{
    bw2=bw;
    bw2(1:2:end,1:2:end)=0;
    [row, col]=find(bw2);
    %}
    % bw2=redundant(bw);
    % [row, col]=find(bw2);
    % 

    %
    for mi=1:d1:m-s1+1
        for ni=1:d2:n-s2+1
            if(mi+s1-1>m||ni+s2-1>n)
                continue;
            end
            mii=(mi-1)*mag+1;
            nii=(ni-1)*mag+1;
            k=k+1;
            tmpl=reshape(img(mi:mi+(s1-1),ni:ni+(s2-1)),s1*s2,1);
            tmph=reshape(img2(mii:mii+(ss1-1),nii:nii+(ss2-1)),ss1*ss2,1);
%             mvalue=mean(tmpl);
%             tmpl=tmpl-mvalue;
%             tmph=tmph-mvalue;
            patch1(:,k)=tmpl;  % reshape(img(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);   % [reshape(img(mi:mi+s1-1,ni:ni+s2-1),s1*s2,1);mi;ni];
            patch2(:,k)=tmph; % reshape(img2(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1); % [reshape(img2(mi:mi+s1-1,ni:ni+s2-1),s1*s2,1);mi;ni];


        end
    end
    %{
    for ii=1:length(row)
            mi=row(ii);
            ni=col(ii);
            if(((mi-(s1-1)/2)<1)||((mi+(s1-1)/2)>m)||((ni-(s2-1)/2)<1)||((ni+(s2-1)/2)>n))
                continue;
            end
            temp=bw(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2);
            temp=double(temp);
            if((sum(sum(temp)))>=3) % s1=5: 0.338 s1=7: 0.242 % s1=10: 0.166 % s1=9: 0.1886
            k=k+1;
            tmpl=reshape(img(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            tmph=reshape(img2(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            mvalue=mean(tmpl);
            tmpl=tmpl-mvalue;
            tmph=tmph-mvalue;
            patch(:,k)=tmpl;  % reshape(img(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);   % [reshape(img(mi:mi+s1-1,ni:ni+s2-1),s1*s2,1);mi;ni];
            patch2(:,k)=tmph; % reshape(img2(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1); % [reshape(img2(mi:mi+s1-1,ni:ni+s2-1),s1*s2,1);mi;ni];
            tmp1=reshape(imgf1(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            tmp2=reshape(imgf2(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            tmp3=reshape(imgf3(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            tmp4=reshape(imgf4(mi-(s1-1)/2:mi+(s1-1)/2,ni-(s2-1)/2:ni+(s2-1)/2),s1*s2,1);
            imgfeature(:,k)=[tmp1;tmp2;tmp3;tmp4];
            % Normalize the patch, patch2 
            %
            % dh=max(patch(:,k));
            % dl=min(patch(:,k));
            % ddh=max(patch2(:,k));
            % ddl=min(patch2(:,k));
            % rate1=255/(dh-dl);
            % rate2=255/(ddh-ddl);
            % patch(:,k)=(patch(:,k)-patch(25,k));
            % patch2(:,k)=(patch2(:,k)-patch2(25,k));
            % [tmp,dd,var]=normalization(patch(:,k));
%             var=norm(patch(:,k),1)/(s1*s2); %sqrt(sum(patch(:,k).^2)/s1/s2);
%             patch(:,k)=patch(:,k)/var;
            % patch2(:,k)=patch2(:,k)-mean(patch2(:,k));
%             patch2(:,k)=patch2(:,k)/var;
            % [patch2(:,k),dd2,var2]=normalization(patch2(:,k));
            
            end
    end
    %}
    % All pacthes in different images are saved in the same variable name
    % 'patch', pls be carefull
%     save([output,flist(i).name(1:end-4)],'patch1');
%     save([output2,flist(i).name(1:end-4)],'patch2');

     L=[L,patch1];
     H=[H,patch2];
    fprintf('N=%d\n',i);
end
