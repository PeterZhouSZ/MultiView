clear all;
str_Load= 'F:\HZ\neighbor embedding\Face_align\database\Facial Expression\Interpolated\5by5\'; %Patch_Learning\Adjust  % 'I:\TraingSet\CMU_Database\Testing_1\Chai\Left+_Left-\20_15\'
str_Load_1='F:\HZ\neighbor embedding\Face_align\database\Facial Expression\phase_3_adjustment\100by100\\'; %30_30\Adjust\ %10_10\Adjust\  Our\Optical-adjust-250  Our\Optical_adjust_t\
for pp = 1:68
str_ref=strcat(str_Load, num2str(pp),'.bmp');    
str_test=strcat(str_Load_1, num2str(pp),'.bmp'); 
I=imread(str_test);
error(pp)=calculate_error(I,'MSE',str_ref);
error3(pp)=calculate_error(I,'MSSIM',str_ref);
end
10*log10(255*255/mean(error))
mean(error3)


error=error';
error3=error3';
x=[1:1:67];
x1=[1:1:68];
figure(1)
plot(x1,error,'r');hold on;
figure(2)
plot(x1,error,'r');hold on;

error_total=[error_total;error,0,error2,0,error3];
save([str_Load_1,'error_total_TIP'],'error_total');

str_Load='E:\TraingSet\CMU_Database\pie_jpg\nov_2000-dec_2000\040';
str_out='D:\HZ\neighbor embedding\Right_22.5\';
for i=55:69
%     if i==33 || i==38
%         continue
%     else
 str_test=strcat(str_Load, num2str(i),'\expression\N_W_05.jpg');
 str_ref=strcat(str_out, num2str(i),'.bmp');
 I=imread(str_test);
 I1=rgb2gray(I);
 imwrite(I1,str_ref);
%     end
end