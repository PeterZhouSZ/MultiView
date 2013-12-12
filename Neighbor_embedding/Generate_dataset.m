clear all;
input='D:\HZ\neighbor embedding\Frontal\Patch\'; %'C:\Documents and Settings\Administrator\My Documents\MATLAB\two_step_SR\H-patch\';
input2='D:\HZ\neighbor embedding\Left_22.5\Patch\'; %'C:\Documents and Settings\Administrator\My Documents\MATLAB\two_step_SR\L-patch\';
output='D:\HZ\neighbor embedding\Frontal\Patch\'; %'C:\Documents and Settings\Administrator\My Documents\MATLAB\two_step_SR\H\';
output2='D:\HZ\neighbor embedding\Left_22.5\Patch\'; %'C:\Documents and Settings\Administrator\My Documents\MATLAB\two_step_SR\L\';
% output3='D:\Huyu\experiment3\Classified\';
extension='mat';
% Generate the whole patch column for high-res.
flist=dir([input,'*.',extension]);
[mm,nn]=size(flist);
H=[];
H=single(H);
L=[];
FL=[];
L=single(L);
FL=single(FL);

for pnum=1:mm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([input,flist(pnum).name]);
    H=[H,patch2];
%     load([input2,flist(pnum).name]);
%     L=[L,patch];
end
save([output,'Hseries'],'H')
% save([output2,'Lseries_z'],'L')
% clear H;
% Generate the whole patch column matrix for low-res.

L=[];
FL=[];
L=single(L);
FL=single(FL);


for k=1:mm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load feature for ¡°Super-resolution through neighbor embedding
   load([input2,flist(k).name]);
    L=[L,patch1];
%     FL=[FL,imgfeature]; % for neighbor embedding
    
          % L=[L,patchL];
end
% clear L2examplar
save([output2,'Lseries'],'L')

% clear L;
% Generate the similarity matrix
% Compute each patch group distance matrix
% S
%{
    k=size(H,2);
    S=zeros(k,k,'single');
    for i=1:k
       for j=1:k
          if(i==j)
              S(i,j)=0;
              continue;
          end
          S(i,j)= -(norm(H(:,i)-H(:,j)))^2;%abs(H(:,i)'*H(:,j))/norm(H(:,i))/norm(H(:,j))-1;  
       end
    end
% Compute the preference
    r=matmedian(S);
    for i=1:k
        S(i,i)=r;
    end
% Affinity propagation
    p=diag(S);
    [idx,netsim,expref,dpsim] = apclustersparse(S,p,'maxits',1000,'dampfact',0.8,'nonoise');
    clear S;
% Classify the samples into several classes
% FOR L
[m,n]=size(L);
i=[1:n];
ind=find(~(double(idx(i,end))-i'));
% EXAMPLARS
k=length(ind);
Lexamplar=L(:,ind);
Lsample=cell(2,k);
for i=1:k
    index=find(~(idx1-ind(i)));
    Lsample(1,i)={L(:,index')};
    Lsample(2,i)={index};
end
save([output3,'LCLASS'],'ind','Lexamplar','Lsample'); 
%}
