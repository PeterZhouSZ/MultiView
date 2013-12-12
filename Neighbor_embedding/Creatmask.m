function temp=Creatmask(img,s1,s2,list)
[m,n]=size(img);
temp=zeros(m,n);
pad=ones(s1,s2);
x=list(1:2:end);
y=list(2:2:end);
K=length(x);
for i=1:K
    temp(x(i):x(i)+s1-1,y(i):y(i)+s2-1)=temp(x(i):x(i)+s1-1,y(i):y(i)+s2-1)+pad;
end
[p,q]=size(temp);
%{
for i=1:p
    for j=1:q
        if(temp(i,j)>0)
            temp(i,j)=temp(i,j)^(-1);
        end
    end
end
%}
ind=temp>0;
temp(ind)=temp(ind).^(-1);
