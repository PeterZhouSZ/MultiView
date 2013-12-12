profile_name = 'X:/MultiView/YaleB_exteneted/left_45/interpolated/4by4/';
frontal_name = 'X:/MultiView/YaleB_exteneted/frontal/phase1/';
pts = [1 14;60 14;30 40; 30 5];
window_size = [9 9];
L =[]; H = [];
for iter = 1:4
    p = pts(iter,:);    
    for frame = 2:11
        imL = imread([profile_name, sprintf('%02d.png', frame)]);
        imH = imread([frontal_name, sprintf('%02d.png', frame)]);
        pL = imcrop(imL, [p window_size]);
        pH = imcrop(imH, [p window_size]);
        L = [L pL(:)];
        H = [H pH(:)];
    end
end

vert1 =  zeros(40, 1, 3);
upper1 = zeros(1, 12, 3);
bottom1 = zeros(1,12,3);
upper1(:,:,3) = 255;
bottom1(:,:,3) = 255;
vert1(:,:,3) = 255;
im = [upper1; vert1 repmat(L, [1 1 3]) vert1; bottom1];
imwrite(im, 'inputPatch.png');

patchTotal = [];
img_size = 10;
width = 1;
upper1 =  zeros(width, 10 * img_size + 2*width + width*(img_size - 1), 3);
upper2 = zeros(width, 10 * img_size + 2*width + width*(img_size - 1), 3); 
vert1 =  zeros(img_size, width, 3);
vert2 =  zeros(img_size, width, 3);
bottom1 = zeros(width, 10 * img_size + 2*width + width*(img_size - 1), 3);
bottom2 = zeros(width, 10 * img_size + 2*width + width*(img_size - 1), 3);

upper1(:,:,1) = 255;
upper2(:,:,2) = 255;
vert1(:,:,1) = 255;
vert2(:,:,2) = 255;
bottom1(:,:,1) = 255;
bottom2(:,:,2) = 255;
demoPatch = [];
for j = 1:4
    patchTotal = [];
    for i = 1:10
        patch1 = L(:,(j - 1)*10 + i);
        patch1 = reshape(patch1, 10 ,10);
        patch1 = repmat(patch1, [1 1 3]);
        patchTotal = [patchTotal vert1 patch1];
    end
    %patchTotal = repmat(patchTotal, [1 1 3]); 
    demoPatch = [demoPatch; patchTotal vert1; bottom1];
end
demoPatch = [upper1; demoPatch; bottom1;];
imshow(uint8(demoPatch))
imwrite(uint8(demoPatch), 'LRpatch.png')

demoPatch2 = [];
for j = 1:4
    patchTotal = [];
    for i = 1:10
        patch1 = H(:,(j - 1)*10 + i);
        patch1 = reshape(patch1, 10 ,10);
        patch1 = repmat(patch1, [1 1 3]);
        patchTotal = [patchTotal vert2 patch1];
    end
    %patchTotal = repmat(patchTotal, [1 1 3]); 
    demoPatch2 = [demoPatch2; patchTotal vert2; bottom2];
end
demoPatch2 = [upper2; demoPatch2; bottom2;];
imshow(uint8(demoPatch2))
imwrite(uint8(demoPatch2), 'HRpatch.png')