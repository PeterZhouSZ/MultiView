% We already have the img, so we want to generate the degradation matrix
function D=genDmat(img)

[m,n]=size(img);
mag=4; % the size of each dimension of img should be equal to mag*integer 
m2=ceil(m/mag);
n2=ceil(n/mag);
% Various filter kernel can be selected
dx=3;
dy=3;
h=fspecial('gaussian',[2*dx+1,2*dy+1],(2*dx+1)/2/3); % Gaussian kernel ...
% h=zeros(2*dx+1,2*dy+1);
% h(dx+1,dy+1)=1;          % Decimation
%
D=sparse(m2*n2,m*n);
p=0;
for j=1:mag:n
    for i=1:mag:m
        tmp=zeros(m,n);
        if(-dx+i<1||dx+i>m||-dy+j<1||dy+j>n)
            tmp(i,j)=1;
        else
            tmp(-dx+i:dx+i,-dy+j:dy+j)=h;
        end
        tmp=tmp(:);
        ind=tmp~=0;
        p=p+1;
        D(p,ind)=tmp(ind);
    end
end
