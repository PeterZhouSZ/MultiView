function imout=recover3(img,outpatch,list,template,s1,s2)
x=list(1:2:end);
y=list(2:2:end);
[m,n]=size(img);
img2=zeros(m,n);
K=length(x);
for i=1:K
    patch=reshape(outpatch(:,i),s1,s2);
    img2(x(i):x(i)+s1-1,y(i):y(i)+s2-1)=img2(x(i):x(i)+s1-1,y(i):y(i)+s2-1)+template(x(i):x(i)+s1-1,y(i):y(i)+s2-1).*patch;
end
% if medpass is not used, and img2 is already the intensity, then use
% recover3
imout=img2;
% To confirm whether ind will be all zeros
%{
ind=template==0;
ind=double(ind);
imout=imout+img.*ind;
hh=fspecial('gaussian',[5 5],5/3/2);
hh2=fspecial('gaussian',[5 5],5/3/2);
ind=imfilter(ind,hh,'replicate');
tmp=(ind>0).*(ind<1);
ind=logical(tmp);
tmp=imout;
tmp=imfilter(tmp,hh2,'replicate');
imout(ind)=tmp(ind);
%}