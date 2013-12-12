% Illumination normalization Xie's method
% img is the input image, N is the size of the local window for each dimension
% LM is the map of local mean, LV is the map of local variance
function [LM,LV]=local_mean_var(img,N)
[m,n]=size(img);
LM=ones(m,n);
LV=ones(m,n);
nn=(N-1)/2;
imgp=padarray(img,[nn,nn],'replicate','both');
for i=nn+1:m+nn
    for j=nn+1:n+nn
        tmp=imgp(i-nn:i+nn,j-nn:j+nn);
        lm=mean(tmp(:));
        lv=std(tmp(:))+0.5;
        LM(i-nn,j-nn)=lm;
        LV(i-nn,j-nn)=lv;
    end
end
        