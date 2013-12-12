%% demo eigentransformation
inDirL = 'X:/MultiView/YaleB_exteneted/right_45/downSample/4by4/';
inDirH = 'X:/MultiView/YaleB_exteneted/right_45/phase1/';
outdir = 'X:\MultiView\YaleB_exteneted\right_45\Results\Tensor\4by4\';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
NumOfSamples = 27;
for i = 1:NumOfSamples
    eigTrans(i,inDirL,inDirH,outdir,NumOfSamples);
end
