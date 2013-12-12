function eigTrans(INDEX,inDirL,inDirH,outdir,NumOfSamples)
% inDirL = '';
% inDirH = '';
% outdir = '';


for i = 1: NumOfSamples

    subjID = i;
    extName = sprintf('%02d.png',subjID); %strcat(str_H_Path,num2str(i), '.bmp');   
    strL = fullfile(inDirL,extName);
    strH = fullfile(inDirH,extName);
    Image = imread(strH);
    Image=double(Image);
    [rw cl] = size(Image);
    Image1= imread(strL);
    Image1=double(Image1);

    TrainingImage(:,i)=double(Image(:));
    TrainingImage_L(:,i)=double(Image1(:));

end
%Low resolution Image

    testFace = TrainingImage_L(:,INDEX);
    TrainingImage_L(:,INDEX) = [];
    TrainingImage(:,INDEX) = [];


    MeanFace = mean(TrainingImage');
    MeanFace_L = mean(TrainingImage_L');

    DemeanFace_L = TrainingImage_L - repmat(MeanFace_L',1, NumOfSamples - 1);
    DemeanFace = TrainingImage - repmat(MeanFace',1, NumOfSamples - 1);

    LOW_A = DemeanFace_L'*DemeanFace_L;
    [V_L D_L] = eig(LOW_A);
    D_L = diag(D_L);

    for i = 0 : 15
        V_L1(:,i+1) = V_L(:,NumOfSamples-1-i);
        D_L1(i+1,1) = D_L(NumOfSamples-1-i,1);
    end

    D_L1 = D_L1.^-0.5;
    D_L1 = diag(D_L1);
    EV_L = DemeanFace_L*V_L1*D_L1;
    %%%Testing


    testFace = testFace - MeanFace_L';
    coef_L = EV_L'*testFace;
    C_L = V_L1*D_L1*coef_L;
    reFace = DemeanFace*C_L + MeanFace';
    reFace = reshape(reFace,rw,cl);

    extName = sprintf('%02d.png',INDEX);
    svName = fullfile(outdir,extName);
    imwrite(uint8(reFace),svName);
