function W=MyLLR2(x,y,B)
% x is D¡Á1 input data, y is D¡ÁK reference data. B is diagonal matrix whoes entries are the weights for each of D dimensions 
if(nargin>2)
% B=eye(size(B));
A=diag(B);
A=A.^0.5;
A=diag(A);
x=x(:);
x=A*x;
y=A*y;
end

% D=size(y,3);
% m=size(y,1);
% n=size(y,2);
% K=m*n;
% xx=1:n;
% yy=1:m;
% [xx,yy]=meshgrid(xx,yy);
% j=0;
% u=zeros(D,K-1);
% k=0;
% for i=1:m
%     for j=1:n
%         if(i==2&&j==2)
%             continue;
%         end
%         k=k+1;
%         u(:,k)=y(i,j,:);
%     end
% end

u=y;

%{
for i=1:K
    if(i==(K+1)/2)
        continue;
    end
    j=j+1;
    u(:,j)=y(yy(i),xx(i),:);
end
%}
%{
C=u'*u;
%
ww=zeros(m,n);
if(K-1>D)
    disp('SINGULAR')
    C=C+0.001*trace(C)*eye(K-1);
end
%
Ci=inv(C);
% alpha=1-sum(sum(Ci.*repmat((x(:)'*u),[K-1,1]),1),2);
alpha=1-sum(Ci*u'*x);
beta=sum(sum(Ci));
lambda=alpha/beta;
w=Ci*(u'*x+lambda);
%}
%

K=size(u,2);
% Use LLE to reconstruct the input with respect to samples
%
z = u-repmat(x,1,K); % shift ith pt to origin
C = z'*z;                                        % local covariance

% ww=zeros(m,n);

% C=single(C);
if(sum(abs(C(:)))==0)
    % disp('singular');
    tol=1e-3;
    C = C + eye(K,K)*tol*1;
end
% 
while(abs(cond(C,1)>200))
%      disp('singular');
     tol=1e-3;
     C = C + eye(K,K)*tol*trace(C);
end

% if(abs(cond(C)>10))
%      disp('singular');
%      tol=1e-3;
%      C = C + eye(K,K)*tol*trace(C);
% end

% C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
W = C\ones(K,1);                           % solve Cw=1
% W=pinv(C)*ones(K,1);
W = (W)/(sum(W)+eps);                  % enforce sum(w)=1
%}

% Use the sparse regression to reconstruct the input
%{
lam=0.5;
cvx_begin
    variable W(K,1);
    minimize( norm(x-u*W,2)+lam*norm(W,1) );
cvx_end;
W=W;
%}    

%}
%{
if(abs(sum(W))-1<0.0001)
    ww=ww;
else
    ww(2,2)=1;
    return;
end
%}
%}

%{
j=0;
for i=1:K
    if(i==(K+1)/2)
        continue;
    end
    j=j+1;
    ww(yy(i),xx(i))=w(j);
end
%}

% k=0;
% for i=1:m
%     for j=1:n
%         if(i==2&&j==2)
%             continue;
%         end
%         k=k+1;
%         ww(i,j)=W(k);
%     end
% end