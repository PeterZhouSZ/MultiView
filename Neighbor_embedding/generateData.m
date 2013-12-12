function [L, H, FL] = generateData(inDirL, inDirH, s1, s2, overlap1, overlap2, mag)

extension='.png';
flist=dir([inDirL,'*',extension]);

ss1 = s1 * mag;
ss2 = s2 * mag;
d1 = s1 - overlap1;
d2 = s2 - overlap2;

h1=[-1 0 1];
h2=h1';
hh1=[1 0 -2 0 1];
hh2=hh1';

L = [];
H = [];
FL = [];

for ii = 1:length(flist)
    k = 0;
    imgL = double(imread([inDirL sprintf('%02d',ii) extension]));
    imgH = double(imread([inDirH sprintf('%02d',ii) extension]));

    [m, n] = size(imgL);
    [mH, nH] = size(imgH);
    
    imgf1=imfilter(imgL,h1,'replicate');
    imgf2=imfilter(imgL,h2,'replicate');
    imgf3=imfilter(imgL,hh1,'replicate');
    imgf4=imfilter(imgL,hh2,'replicate');
    
    for mi=1:d1:m-s1+1
        for ni=1:d2:n-s2+1
            mii=(mi-1)*mag+1;
            nii=(ni-1)*mag+1;
            if(mi+s1-1>m || ni+s2-1>n || mii+ss1-1 > mH ||nii+ss2-1 > nH)
               continue;
            end
            %{
            temp=bw(mi:mi+s1-1,ni:ni+s2-1);
            temp=double(temp);
            h = fspecial('gaussian',[s1 s1],s1/2/3);
            %}

            k=k+1;
            patch1(:,k) = reshape(imgL(mi:mi+s1-1,ni:ni+s2-1),(s1)*(s2),1);
            pm(k) = mean(patch1(:,k));
            patch1(:,k) = patch1(:,k) - pm(k); 
            
            patch2(:,k) = reshape(imgH(mii:mii+ss1-1,nii:nii+ss2-1),(ss1)*(ss2),1);
            
            tmp1=reshape(imgf1(mi:mi+s1-1,ni:ni+s2-1),(s1)*(s2),1);
            tmp2=reshape(imgf2(mi:mi+s1-1,ni:ni+s2-1),(s1)*(s2),1);
            tmp3=reshape(imgf3(mi:mi+s1-1,ni:ni+s2-1),(s1)*(s2),1);
            tmp4=reshape(imgf4(mi:mi+s1-1,ni:ni+s2-1),(s1)*(s2),1);
            
            imgfeature(:,k)=[tmp1;tmp2;tmp3;tmp4];
            
            
            % Contrast normalization
            % dl(k)=min(inputpatch(:,k));
            % dh(k)=max(inputpatch(:,k));
            % rate(k)=255/(dh(k)-dl(k));
            % dmean(k)=inputpatch(25,k);
            % inputpatch(:,k)=inputpatch(:,k)-dmean(k);
            % [tmp,dd(k),var(k)]=normalization(inputpatch(:,k));
    %         pm(k)=mean(inputpatch(:,k));
    %         var(k)=std(inputpatch(:,k));                     % sqrt(sum(inputpatch(:,k).^2)/s1/s2);norm(inputpatch(:,k),1)/((s1+2)*(s2+2));
    %         inputpatch(:,k)=inputpatch(:,k)-pm(k);

    % inputpatch(:,k)=inputpatch(:,k)/var(k);

        end
    end
    L = [L patch1];
    H = [H patch2];
    FL = [FL imgfeature];
end