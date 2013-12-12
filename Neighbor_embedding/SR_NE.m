% function inputdeal(img,L1,L2,L3,L4)
% Reference: H. Chang, D.-Y. Yeung, and Y. Xiong, ¡°Super-resolution through neighbor embedding,¡± in Proc. IEEE Conf. Computer Vision and Pattern Recognition, 2004, pp. I-275-I-282.
%% All the rights reserved by zhuo hui, Carnegie Mellon University
clear;
% Hpath, Lpath are used in line 127,128. They are used to load the
% generated database (LR-HR patch pairs). The generation of database is described in the reference
Lpath = 'X:/MultiView/YaleB_exteneted/right_45/downSample/4by4/';
Hpath=  'X:/MultiView/YaleB_exteneted/right_45/phase1/';
% FLpath='E:\HZ\neighbor embedding\Frontal View_Novel\Eigen_Frontal\Eigen_4\140_209\Feature\';
imout_path='X:\MultiView\YaleB_exteneted\right_45\Results\NE\4by4\';% path for output

if ~exist(imout_path, 'dir')
    mkdir(imout_path);
end
% img_path='E:\HZ\neighbor embedding\Frontal View_Novel\Eigen_Frontal\Eigen_4\140_209\'; % path for LR image
% % img_path_hr='F:\TraingSet\CMU_Database\Testing\Frontal\';% path for HR image which is used to evaluate the reconstruction error
extension='.png';
flist=dir([Hpath,'*',extension]);
% [pp,qq]=size(flist);
% list=cell(pp,1);
% outpatch=cell(pp,1)
mag = 4;

s1 = 7;   %7
s2 = 7;
ss1 = s1*mag;
ss2 = s2*mag;

overlap1 = 6;   %6
overlap2 = 6;
d1 = s1-overlap1;
d2 = s2-overlap2;


% h1 = [-1 0 1];
% h2 = h1';
% hh1 = [1 0 -2 0 1];
% hh2 = hh1';
[L, H, feat] = generateData(Lpath, Hpath, s1, s2, overlap1, overlap2, mag);
numOfpatch = size(H,2)/length(flist);
inPatch = [];

outPatch=[];
inputfeature=[];
% Batch processing for a series of test photo

for pnum = 1 : length(flist)
    img_name = [Lpath sprintf('%02d',pnum) extension];%strcat(Lpath,num2str(pnum),extension);
    img=imread(img_name);
    [m, n] = size(img);
    img=double(img);
    imtmp=img;
   
    L1 = L;
    H1 = H;
    featL = feat;
    
    strPatch = (pnum-1) * numOfpatch + 1;
    endPatch = pnum * numOfpatch;
    inPatch = L1(:,strPatch : endPatch); 
    inputfeature = featL(:,strPatch : endPatch);
    L1(:,strPatch : endPatch) = [];
    H1(:,strPatch : endPatch) = [];
    featL(:,strPatch : endPatch) = [];
    list=[];
   for mi=1:d1:m-s1+1
    for ni=1:d2:n-s2+1
        if(mi+s1-1>m||ni+s2-1>n)
           continue;
        end
        mii=(mi-1)*mag+1;
        nii=(ni-1)*mag+1;
        list=[list,[mii,nii]];
    end
   end

% Compute input's feature

   
% Normalize the inputpatch
% inputpatch=inputpatch/255;

% Process inputpatch
% load([Lgpath,'L1']);
% [m,n]=size(inputpatch);
% k=6; % This should be consistent with the K in the operate2 function
% w=zeros(k,n);
% Label=operate2(input_L1,ind1,L1examplar,L1sample);
% for i=1:n
%     P=L1(:,Label(i,:));
%     w(:,i)=P\input_L1(:,i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Search the best matching pathes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% Generate adj
    adj=adjacent(inPatch,list,ss1,ss2);
    % There are two optional selection below:
    % outpatch=search(inputpatch,H,ind,Lexamplar,Lsample);
    K = 5;
% [Label,B]=Knn(inputpatch,L,K);


    [Label,B]=knnsearch(inputfeature',featL',K);
    B=B.*B;
% 
    for j=1:size(inputfeature,2)
        HK=H1(:,Label(j,:));
        LK=L1(:,Label(j,:));
        ww = MyLLR2(inPatch(:,j),LK);
        outPatch(:,j) = HK * ww;
        outPatch(:,j) = outPatch(:,j);
    end
%     R = recover(outPatch,ss1,ss2,100,100,overlap1);
%     imshow(R, [])

% clear L;
% adj=logical(adj);
% clear img2;
% clear img;
% clear inputpatch;
% % Generate lambda matrix 
% lambda=Glambda(H,Label,adj,list,K,s1,s2);
% 
% % Generate local matrix
% local=cell(1,size(adj,2));
% [bm,bn]=size(B);
% T=8000;
% for iii=1:bm
%     local{iii}=B(iii,:)';
%     local{iii}=exp(-1*local{iii}/T);
% end
% % clear B;
% % belief propagation
% [bel, convereged, I_new]=loopybp(adj,lambda,local,K);
% % Generate the outpatch
% for jjj=1:bm
%     outpatch(:,jjj)=H(:,Label(jjj,I_new(jjj)));
% end


%Label=operate(inputpatch,ind,Lexamplar,Lsample);
%[mp,np]=size(Label);
%for i=1:mp
%    P=L(:,Label(i,:));
%    w(:,i)=P\inputpatch(:,i);
%end

%for i=1:mp
%    P=H(:,Label(i,:));
%    outpatch=[outpatch,P*w(:,i)];
%end

% Denormalize the outpatch
% outpatch=outpatch*255;
% for i=1:length(var)
%     outpatch(:,i)=outpatch(:,i)*var(i);
% end        
% for i=1:length(pm)
%     outpatch(:,i)=outpatch(:,i)+pm(i);
% end

%%%%%%%%%%%
% Recovery
%%%%%%%%%%%
    imtmp=imresize(imtmp,size(imtmp)*mag);
    template=Creatmask(imtmp,ss1,ss2,list);
    imout=recover3(imtmp,outPatch,list,template,ss1,ss2);
    
%     imout=imresize(imout,[124 108],'bilinear');
%       str_out=strcat(imout_path,num2str(pnum),extension);
%       str_hr=strcat(img_path_hr,num2str(pnum),extension);
%     imwrite(imout,str_out);
%     error(pnum)=calculate_error(double(imout),'MSE',str_hr);
%     error2(pnum)=calculate_error(double(imout),'ESMSE',str_hr);
%     error3(pnum)=calculate_error(double(imout),'MSSIM',str_hr);
    imout=uint8(imout);
    str_out=strcat(imout_path,num2str(pnum),extension);
    imwrite(imout,str_out);
    fprintf('N=%d\n',pnum);
%     imshow(imout,[0 255]);
    
end
% save([imout_path,'error_total'],'error','error2','error3');