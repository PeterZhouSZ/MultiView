% Face
% Functions used
%   NormalizeToEyes_T()
%   NormalizeImageByMeanVar()
% by Thomas
% Last updated on 30 Jun 2007
function faceNormalize(filepath_image, output1, width, height, r1, r2, r3)
% The parameters are used to adapt the ratio of the height and width of 
% normalized facial images

%r1 = 1.35;  %ORT 1.40  % 1.8 %2 scaling image %1.45  %1.6 left_22.5   %left_45 1.57  %1.6
%r2 = 0.22;  %ORI 0.23  0.26 %up and down up: -   %0.25  
%r3 = 0.09;  %ORI 0.10 %0.2  left_right expand -    %0.15  %0.11 left_22.5

%width = 100;
%height = 100;
nC = 68;% no. of classes, max. 169 %2

% read the images
% filepath_person_template = 'I:\\FERET_CD1\\partitions\\by_subject\\%05d_basenames.txt';
% filepath_eye_template = 'D:\\Huyu\\database\\GTdb_crop\\labels_gt\\labels\\lab%03d.txt';

%filepath_image = 'E:\HZ\neighbor embedding\Face_align\database\phase_2_adjustment\';        %  path for the read-in facial images
%output1='E:\HZ\neighbor embedding\Face_align\database\phase_3_adjustment\';         %  path for the output normalized facial images
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

list=dirlist.list;
for ii = 1 : length(list)   %%%%%%%%%%%
    eye_pos1=[list(ii).Leyex,list(ii).Leyey];
    eye_pos2=[list(ii).Reyex,list(ii).Reyey];
    eye_pos2=double(eye_pos2);
    eye_pos1=double(eye_pos1);
    img=imread([filepath_image,list(ii).name]);
    if(size(img,3)>1)
        img=rgb2gray(img);
    end
    img=double(img);
%     [height, width]=size(img);
    B = NormalizeToEyes_T(img,eye_pos1(1),eye_pos1(2),eye_pos2(1),eye_pos2(2),r1,r2,r3,width,height);
    imwrite(uint8(B),[output1,list(ii).name]);
    fprintf('N = %d\n',ii)% for display
end
