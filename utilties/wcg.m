% pixel structures learning 
function [B,C]=wcg(A,sigma,s1,b,ermap)
% % Input: A: is a matrix of size m*n*K, which stores the reference examples. (m,n) is the size of the image, K is the number of reference examples. 
%          sigma: is a parameter, see our paper. 
%          s1,b:  are used in case of learning the local pixel structures of a patch (in the future version). Now, we simply set s1=2, b=1. 
%          ermap: is the matrix that stores the warping error maps. Its size: m*n*K
% Output:  B: the learned embedding weights. B{1,1} constains all the embedding weights corresponding to the left-up direction connection between two neighboring pixels...
%          C: the corresponding proportion factor. 
% 

sb=s1-b;
eps=1;
% 8 directions
[m,n,p]=size(A);
B=cell(3,3);
a=2;
sigma=sigma^a;

%{
B{1,1}=exp(-sum((A(2:end-1,2:end-1,:)-A(1:end-2,1:end-2,:)).^a/sigma/size(A,3),3));
B{1,2}=exp(-sum((A(2:end-1,2:end-1,:)-A(1:end-2,2:end-1,:)).^a/sigma/size(A,3),3));
B{1,3}=exp(-sum((A(2:end-1,2:end-1,:)-A(1:end-2,3:end,:)).^a/sigma/size(A,3),3));
B{2,1}=exp(-sum((A(2:end-1,2:end-1,:)-A(2:end-1,1:end-2,:)).^a/sigma/size(A,3),3));
B{2,2}=0*exp(-sum((A(2:end-1,2:end-1,:)-A(2:end-1,2:end-1,:)).^a/sigma/size(A,3),3)); %zeros(m-2,n-2);  % ones(m-2,n-2); 
B{2,3}=exp(-sum((A(2:end-1,2:end-1,:)-A(2:end-1,3:end,:)).^a/sigma/size(A,3),3));
B{3,1}=exp(-sum((A(2:end-1,2:end-1,:)-A(3:end,1:end-2,:)).^a/sigma/size(A,3),3));
B{3,2}=exp(-sum((A(2:end-1,2:end-1,:)-A(3:end,2:end-1,:)).^a/sigma/size(A,3),3));
B{3,3}=exp(-sum((A(2:end-1,2:end-1,:)-A(3:end,3:end,:)).^a/sigma/size(A,3),3));
%}

% For LLE_try_block
%{
weight=1./(ermap+eps).^2;
weight=reshape(weight,[1,1,p]);
fen=sum(weight);
weight=repmat(weight,[m-2*sb,n-2*sb,1]);
%}
% For the global image
%
ddx=9;
ddy=9;
h=ones(ddx,ddy)/ddx/ddy;
for jj=1:size(ermap,3)
    ermap(:,:,jj)=imfilter(ermap(:,:,jj),h,'replicate');
end
ermap=ermap(1+sb:end-sb,1+sb:end-sb,:);
weight=1./(ermap+eps).^2;
fen=sum(weight,3);
fen=repmat(fen,[1,1,size(weight,3)]);
%}

%{
B{1,1}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb:end-sb-1,:)).^a/sigma/size(A,3),3));
B{1,2}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb+1:end-sb,:)).^a/sigma/size(A,3),3));
B{1,3}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb+2:end-sb+1,:)).^a/sigma/size(A,3),3));
B{2,1}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,sb:end-sb-1,:)).^a/sigma/size(A,3),3));
B{2,2}=0*exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,1+sb:end-sb,:)).^a/sigma/size(A,3),3)); %zeros(m-2,n-2);  % ones(m-2,n-2); 
B{2,3}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,2+sb:end-sb+1,:)).^a/sigma/size(A,3),3));
B{3,1}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,sb:end-sb-1,:)).^a/sigma/size(A,3),3));
B{3,2}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,sb+1:end-sb,:)).^a/sigma/size(A,3),3));
B{3,3}=exp(-sum((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,2+sb:end-sb+1,:)).^a/sigma/size(A,3),3));
%}
B{1,1}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb:end-sb-1,:)).^a).*weight/sigma./fen,3));
B{1,2}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb+1:end-sb,:)).^a).*weight/sigma./fen,3));
B{1,3}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(sb:end-sb-1,sb+2:end-sb+1,:)).^a).*weight/sigma./fen,3));
B{2,1}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,sb:end-sb-1,:)).^a).*weight/sigma./fen,3));
B{2,2}=0*exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,1+sb:end-sb,:)).^a).*weight/sigma./fen,3)); %zeros(m-2,n-2);  % ones(m-2,n-2); 
B{2,3}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(1+sb:end-sb,2+sb:end-sb+1,:)).^a).*weight/sigma./fen,3));
B{3,1}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,sb:end-sb-1,:)).^a).*weight/sigma./fen,3));
B{3,2}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,sb+1:end-sb,:)).^a).*weight/sigma./fen,3));
B{3,3}=exp(-sum(((A(1+sb:end-sb,1+sb:end-sb,:)-A(2+sb:end-sb+1,2+sb:end-sb+1,:)).^a).*weight/sigma./fen,3));
%{
v2=(A(1:end-2,1:end-2,:).*repmat(B{1,1},[1 1 p])+A(1:end-2,2:end-1).*repmat(B{1,2},[1 1 p])+A(1:end-2,3:end).*repmat(B{1,3},[1 1 p])+A(2:end-1,1:end-2).*repmat(B{2,1},[1 1 p])+A(2:end-1,2:end-1).*repmat(B{2,2},[1 1 p])+A(2:end-1,3:end).*repmat(B{2,3},[1 1 p])+A(3:end,1:end-2).*repmat(B{3,1},[1 1 p])+A(3:end,2:end-1).*repmat(B{3,2},[1 1 p])+A(3:end,3:end).*repmat(B{3,3},[1 1 p])+eps);
v1=A(2:end-1,2:end-1,:);
%}
v2=(A(sb:end-sb-1,sb:end-sb-1,:).*repmat(B{1,1},[1 1 p])+A(sb:end-sb-1,sb+1:end-sb,:).*repmat(B{1,2},[1 1 p])+A(sb:end-sb-1,sb+2:end-sb+1,:).*repmat(B{1,3},[1 1 p])+A(sb+1:end-sb,sb:end-sb-1,:).*repmat(B{2,1},[1 1 p])+A(sb+1:end-sb,sb+1:end-sb,:).*repmat(B{2,2},[1 1 p])+A(sb+1:end-sb,sb+2:end-sb+1,:).*repmat(B{2,3},[1 1 p])+A(sb+2:end-sb+1,sb:end-sb-1,:).*repmat(B{3,1},[1 1 p])+A(sb+2:end-sb+1,sb+1:end-sb,:).*repmat(B{3,2},[1 1 p])+A(sb+2:end-sb+1,sb+2:end-sb+1,:).*repmat(B{3,3},[1 1 p])+eps);
v1=A(sb+1:end-sb,sb+1:end-sb,:);
% v3=B{1,1}+B{1,2}+B{1,3}+B{2,1}+B{2,2}+B{2,3}+B{3,1}+B{3,2}+B{3,3};
% C2=1.5./v3;
% C=sum(v2.*v1,3)./(sum(v2.*v2,3)+eps);
C=sum(v2.*(weight./fen).*v1,3)./(sum(v2.*(weight./fen).*v2,3)+eps);
% ind=C>1;
% C(ind)=1;
% ind=C-C2>0;
% C(ind)=C2(ind);
% ind=C>20;
% C(ind)=1;
% pad=ones(m-2,n-2);
% B{2,2}(ind)=pad(ind);
B{1,1}=B{1,1}.*C;
B{1,2}=B{1,2}.*C;
B{1,3}=B{1,3}.*C;
B{2,1}=B{2,1}.*C;
B{2,2}=B{2,2}.*C;
B{2,3}=B{2,3}.*C;
B{3,1}=B{3,1}.*C;
B{3,2}=B{3,2}.*C;
B{3,3}=B{3,3}.*C;