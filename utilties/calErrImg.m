function errImg = calErrImg(strIn, strOut)
recImg = imread(strIn);
Img = imread(strOut);
errImg = abs(recImg - Img);