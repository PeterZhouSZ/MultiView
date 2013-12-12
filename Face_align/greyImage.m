%% convert img format
function greyImage(indir, outdir, indir1, Vind)
% indir = 'E:/TraingSet/CMU_Database/pie_jpg/oct_2000-nov_2000/';
% outdir = 'database/Gray_images/';
index = 1;
for ii = 0:54
    if ii == 33 || ii == 38
        continue;
    else
        filename = sprintf('040%02d/expression/N_W_%02d.jpg',ii, Vind);
        imgname = [indir filename];
        if(exist(imgname))
            imgRGB = imread(imgname);
        else
             filename = sprintf('040%02d/expression/N_G_%02d.jpg',ii, Vind);
             imgname = [indir filename];
             imgRGB = imread(imgname);
        end
        imgG = rgb2gray(imgRGB);
        imwrite(imgG,[outdir num2str(index) '.bmp']);
        index = index + 1;
    end
end

% indir1 = 'E:/TraingSet/CMU_Database/pie_jpg/nov_2000-dec_2000/';
% outdir = 'database/Gray_images/';
index = 54;
for ii = 55:69

        filename = sprintf('040%02d/expression/N_W_%02d.jpg',ii, Vind);
        imgname = [indir1 filename];
        if(exist(imgname))
            imgRGB = imread(imgname);
        else
             filename = sprintf('040%02d/expression/N_G_%02d.jpg',ii, Vind);
             imgname = [indir1 filename];
             imgRGB = imread(imgname);
        end
        imgG = rgb2gray(imgRGB);
        imwrite(imgG,[outdir num2str(index) '.bmp']);
        index = index + 1;

end
