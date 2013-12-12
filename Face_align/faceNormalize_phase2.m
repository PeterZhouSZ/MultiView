function faceNormalize_phase2(filepath_image,output1)

%filepath_image = 'E:\HZ\neighbor embedding\Face_align\database\phase_1_adjustment\';        %  path for the read-in facial images
%output1='E:\HZ\neighbor embedding\Face_align\database\phase_2_adjustment\';         %  path for the output normalized facial images

dirlist=load([filepath_image 'Myannotate.mat']);  % Myannotate.mat contains the positons of eyes for facial images.
% output2='D:\Huyu\database\Feret_CD1\test\';
% filepath_output = 'F:\\temp\\img\\';

% train images
X_train = [];X_test = [];
label_train = [];label_test = [];
N = 1;
% fid = fopen(filepath_train,'r');
% if(fid == -1),fprintf('Error to open the file.\n');end
idx = 1;
filename_train = [];
filename_test = [];
nTrain=0;
nTest=0;
% cl=100;
% rw=100;

list=dirlist.list;

eye_pos_fix_1=[list(1).Leyex,list(1).Leyey];
eye_pos_fix_2=[list(1).Reyex,list(1).Reyey];
mouth_pos_fix=[list(1).Cmouthx,list(1).Cmouthy];

eye_pos_fix_2=double(eye_pos_fix_2);
eye_pos_fix_1=double(eye_pos_fix_1);
mouth_pos_fix=double(mouth_pos_fix);

 x_length_fix=eye_pos_fix_2(1)-eye_pos_fix_1(1);
 y_length_fix=mouth_pos_fix(2)-0.5*(eye_pos_fix_1(2)+eye_pos_fix_2(2));

 y_center_fix=round(0.5*(eye_pos_fix_1(2)+eye_pos_fix_2(2)));

 for ii = 40%:length(list)

    eye_pos1=[list(ii).Leyex,list(ii).Leyey];
    eye_pos2=[list(ii).Reyex,list(ii).Reyey];
    mouth_pos=[list(ii).Cmouthx,list(ii).Cmouthy];
    
    eye_pos2=double(eye_pos2);
    eye_pos1=double(eye_pos1);
    mouth_pos=double(mouth_pos);
    
    img=imread([filepath_image,list(ii).name]);
    img=double(img);
    
    x_length=eye_pos2(1)-eye_pos1(1);
    y_length=mouth_pos(2)-0.5*(eye_pos1(2)+eye_pos2(2));
    
    y_center=round(0.5*(eye_pos1(2)+eye_pos2(2)));
    [rw cl]= size(img);
    
    I1=img(:,eye_pos1(1):eye_pos2(1));
    I2=img(:,eye_pos2(1)+1:end);
    I3=img(:,1:eye_pos1(1)-1);
    I1=imresize(I1,[rw,x_length_fix],'bicubic');
    I2=imresize(I2,[rw,(cl-eye_pos2(1))*x_length_fix/x_length],'bicubic');
    I3=imresize(I3,[rw,(eye_pos1(1)-1)*x_length_fix/x_length],'bicubic');
    [rw,cl3]=size(I3);
    [rw,cl2]=size(I2);
    [rw,cl1]=size(I1);
    img_r(:,1:cl3)=I3;
    img_r(:,cl3+1:cl1+cl3)=I1;
    img_r(:,cl1+cl3+1:cl1+cl2+cl3)=I2;
   
    
    I1_r=img_r(y_center:mouth_pos(2),:);
    I2_r=img_r(mouth_pos(2)+1:end,:);
    I3_r=img_r(1:y_center-1,:);
    I1_r=imresize(I1_r,[y_length_fix,cl1+cl2+cl3],'bicubic');
    I2_r=imresize(I2_r,[(rw-mouth_pos(2))*y_length_fix/y_length,cl1+cl2+cl3],'bicubic');
    I3_r=imresize(I3_r,[(y_center-1)*y_length_fix/y_length,cl1+cl2+cl3],'bicubic');
    
    [rw3,cl_r]=size(I3_r);
    [rw2,cl_r]=size(I2_r);
    [rw1,cl_r]=size(I1_r);
    clear img_n
    img_n(1:rw3,:)=I3_r;
    img_n(rw3+1:rw1+rw3,:)=I1_r;
    img_n(rw1+rw3+1:rw1+rw2+rw3,:)=I2_r;
    
%     new_cl=ceil(cl*x_length_fix/x_length);
%     new_rw=ceil(rw*y_length_fix/y_length);
    
%     I1=img(y_center:mouth_pos(2),:);
%     I1=imresize(I1,[y_length_fix,new_cl],'bicubic');
%     img(y_center:y_center+y_length_fix-1,:)=I1;
    
%      img=imresize(img,[new_rw,new_cl],'bicubic');
%     B = NormalizeToEyes_T(img,eye_pos1(1),eye_pos1(2),eye_pos2(1),eye_pos2(2),r1,r2,r3,width,height);
%     img_n=imcrop(img_n,[rw3-10,cl3-20,100,100]);
    imwrite(uint8(img_n),[output1,list(ii).name]);
    fprintf('N = %d\n',ii)% for display
    clear img_r img_n I1 I2 I3 I1_r I2_r I3_r
end