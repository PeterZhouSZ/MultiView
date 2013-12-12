function downSample(inDir, outDir1, outDir2, factor)
% str_Load='E:\HZ\neighbor embedding\Face_align\database\phase_3_adjustment\';
% str_o='E:\HZ\neighbor embedding\Face_align\database\\LR\';
ext='.png';
for kk = 1 : 27
    %str_H=strcat(inDir,num2str(kk),extension);
     
    str_H = [inDir sprintf('%02d', kk) ext];
    
    str_out1 = [outDir1 sprintf('%02d', kk) ext];
    str_out2 = [outDir2 sprintf('%02d', kk) ext];
    
    Image=imread(str_H);
    [rw, cl]=size(Image);
    dnImage = Image(1:factor:end, 1:factor:end);
    imwrite(uint8(dnImage),str_out1);
    Image = imresize(dnImage, [rw cl], 'bicubic');
    %B_4=genDmat(Image);
    %Image=B_4*(double(Image(:)));    
    %Image=reshape(Image,25,25);
    % imwrite(uint8(Image),str_out);

%     Image=imresize(Image,[100,100],'bicubic');
    imwrite(uint8(Image),str_out2);
end