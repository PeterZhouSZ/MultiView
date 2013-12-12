function outpatch=pca_w(test,imspace,Himspace)
% test=L1(:,index);
% imspace=Neighbors;
% Himspace=Neighbors_H;

psi=mean(imspace,2);
L=imspace-repmat(psi,[1,size(imspace,2)]);
psi_H=mean(Himspace,2);
H=Himspace-repmat(psi_H,[1,size(Himspace,2)]);
test=test-psi;

%
% Constructing the eigenface space. Run this only for the first time. 
[V,D]=eig(L'*L);
diagonal=diag(D);
[tmp,ind]=sort(diagonal,'descend');

dim = 10; %15 %10               % The number of eigenvectors to be used

V=V(:,ind(1:dim));
D=diag(diagonal(ind(1:dim)));
diagonal=diag(D);
%%%%%%%%%%%%%%%%%%%%%%%%% Restrict
% for i=1:dim
%     if diagonal(i)<10
%         INDEX=i-1;
%          break;
%     else
%         continue;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagonal=diagonal.^-0.5;
% if INDEX<0
%   outpatch=Reference;
% else
D2=diag(diagonal);
% V=V(:,ind(1:INDEX));
%{
L=single(L);
V=single(V);
D2=single(D2);
%}
% clear imspace imspace2 D

E=L*V*D2;
w=E'*test;
imspacew=V*D2*w;
outpatch=H*imspacew+psi_H;
% end
