clear all;
load('iter_mse+ssim.mat');
for i=1:11
    iter_mean_1(i)=mean(iter_error_1(i,:));
    iter_mean_2(i)=mean(iter_error_2(i,:));
end
x=[1:10:110];
iter_mean_1=(iter_mean_1-repmat(mean(iter_mean_1),1,length(iter_mean_1)))./(max(iter_mean_1)-min(iter_mean_1)).*5;
iter_mean_2=(iter_mean_2-repmat(mean(iter_mean_2),1,length(iter_mean_2)))./(max(iter_mean_2)-min(iter_mean_2)).*5;
for i=1:20
    iter_error_1(:,i)=(iter_error_1(:,i)-repmat(mean(iter_error_1(:,i)),length(iter_error_1(:,i)),1))./(max(iter_error_1(:,i))-min(iter_error_1(:,i))).*5;
    plot(x,iter_error_1(:,i),'r');
    hold on
end
for i=21:42
    iter_error_1(:,i)=(iter_error_1(:,i)-repmat(mean(iter_error_1(:,i)),length(iter_error_1(:,i)),1))./(max(iter_error_1(:,i))-min(iter_error_1(:,i))).*5;
    plot(x,iter_error_1(:,i),'g');
    hold on
end
for i=43:68
    iter_error_1(:,i)=(iter_error_1(:,i)-repmat(mean(iter_error_1(:,i)),length(iter_error_1(:,i)),1))./(max(iter_error_1(:,i))-min(iter_error_1(:,i))).*5;
    plot(x,iter_error_1(:,i),'b');
    hold on
end
   plot(x,iter_mean_1,'c:','LineWidth',7);   

for i=1:20
    iter_error_2(:,i)=(iter_error_2(:,i)-repmat(mean(iter_error_2(:,i)),length(iter_error_2(:,i)),1))./(max(iter_error_2(:,i))-min(iter_error_2(:,i))).*5;
    figure(2)
    plot(x,iter_error_2(:,i),'r');
    hold on
end
for i=21:42
    iter_error_2(:,i)=(iter_error_2(:,i)-repmat(mean(iter_error_2(:,i)),length(iter_error_2(:,i)),1))./(max(iter_error_2(:,i))-min(iter_error_2(:,i))).*5;
    figure(2)
    plot(x,iter_error_2(:,i),'g');
    hold on
end
for i=43:68
    iter_error_2(:,i)=(iter_error_2(:,i)-repmat(mean(iter_error_2(:,i)),length(iter_error_2(:,i)),1))./(max(iter_error_2(:,i))-min(iter_error_2(:,i))).*5;
    figure(2)
    plot(x,iter_error_2(:,i),'b');
    hold on
end
figure(2)
plot(x,iter_mean_2,'c:','LineWidth',7);

G_T='J:\TraingSet\CMU_Database\Testing\Frontal_4\';
load='J:\TraingSet\CMU_Database\Testing\Frontal_4\iter_reconstruction\sub_';
outdir='J:\TraingSet\CMU_Database\Testing\Frontal_4\iter_reconstruction\error_image';
for i=1:3
    str_ref=strcat(G_T,num2str(i),'.bmp');
    str_load=strcat(load,num2str(i));
    for j=10:10:100
        I=imread(str_ref);
        str_test = strcat(str_load,'_iteration_',num2str(j),'.png');
        I1=imread(str_test);
        error_image = abs(double(I)-double(I1));
        [m,n]=size(error_image);
        for ii = 1:n
        error_image(ii,:) = (error_image(ii,:)-mean(error_image)) ./ (max(error_image)-min(error_image)).*255 ;
        end
        str_iter=strcat('sub_%d_iteration_',num2str(j),'.png');
        savename = sprintf(str_iter,i);
        imwrite(uint8(error_image),fullfile(outdir,savename));
    end
end

    
        
        