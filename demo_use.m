strL= 'E:\HZ\neighbor embedding\Face_align\database\Interpolated\\'; %Patch_Learning\Adjust  % 'I:\TraingSet\CMU_Database\Testing_1\Chai\Left+_Left-\20_15\'
strH='E:\HZ\neighbor embedding\Face_align\database\phase_3_adjustment\\';
strOut = 'E:\HZ\neighbor embedding\func\meanface.png';
for i = 1 : 68
    img = imread([strL sprintf('%d.bmp',i)]);
    imgHR = imread([strH sprintf('%d.bmp',i)]);
    trL(:,i) = img(:);
    trH(:,i) = imgHR(:);
end
meanFL = mean(trL,2);
meanFH = mean(trH,2);

imwrite(uint8(reshape(meanFL,100,100)),strOut)