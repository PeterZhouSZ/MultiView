function Y=recover(outpatch,s1,s2,m,n,overlap)


[m1 n1]=size(outpatch);
m11=1+((n-s1)/(s1-overlap));
% n11=m-s2+1;
Index_Matrix=[];
for v=1:n1
    for u=1:m1
        Index_R=(s1-overlap)*(floor((v-1)/m11))+u-s1*floor((u-1)/s2); %floor((u-1)/s2)+1
        Index_C=(s2-overlap)*(v-1-m11*floor((v-1)/m11))+1+(floor((u-1)/s2));
        Index_Matrix(u,v)= (Index_C-1)*m+Index_R;
    end
end

mark=zeros(m,n);
Recons=zeros(m,n);
mark=reshape(mark,m*n,1);
Recons=reshape(Recons,m*n,1);

for i=1:(m*n)
    for j=1:m1
        for w=1:n1
        if Index_Matrix(j,w)==i
            Recons(i)=Recons(i)+outpatch(j,w);
            mark(i)=mark(i)+1;
           end
        end
    end
end
for i=1:(m*n)
    Y(i)=Recons(i)/mark(i);
end

Y=reshape(Y,m,n);

