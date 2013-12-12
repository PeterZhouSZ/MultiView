function Adjust(im_path,imout_path)
% clear all;
% close all;
% im_path = 'E:\HZ\neighbor embedding\Face_align\database\phase_1_adjustment\';
% imout_path='E:\HZ\neighbor embedding\Face_align\database\phase_1_adjustment\';
extension='png';
flist=dir([im_path,'*.',extension]);
% load U;
[NumOfSample,qq]=size(flist);
x=0;
for pnum=1:NumOfSample
Image_name=flist(pnum).name;  %'00084.bmp';
Image=imread([im_path,Image_name]);
I1=dct2(Image);
x=I1(1,1)+x;
end
x=x/NumOfSample;
for pnum=1:NumOfSample
Image_name=flist(pnum).name;  %'00084.bmp';
Image=imread([im_path,Image_name]);
I1=dct2(Image);
I1(1,1)=x;
I2=idct2(I1);
I2=uint8(I2);
I2=histeq(I2);
I2=imadjust(I2,[],[0.2 0.8],0.8);
imwrite(uint8(I2),[imout_path,Image_name]);
fprintf('N = %d\n',pnum)% for display
end 