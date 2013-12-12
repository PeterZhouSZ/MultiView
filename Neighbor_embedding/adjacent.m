% Notice that this function is used in the new patch grid latice!
function adj=adjacent(inputpatch,list,s1,s2)
% s1=5;
% s2=5;
overlap1=1;
overlap2=1;
d1=s1-overlap1;
d2=s2-overlap2;
[m,n]=size(inputpatch);
% adj=zeros(n,n,'int8');
adj=sparse(n,n);
x=list(1:2:end);
y=list(2:2:end);

for i=1:n
    for j=1:n
        dx=abs(x(i)-x(j));
        dy=abs(y(i)-y(j));
        if(dx<=d1&&dy<=d2)
            adj(i,j)=1;
        end
    end
end
