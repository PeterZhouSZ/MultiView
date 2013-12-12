% img1 is the input image
clear all;
% kk=60;

% KK=9;
sigma2 = 30; %35  Smooth the image the large, the more smoothed results obtained. 
str_i = 'X:/MultiView/YaleB_exteneted/right_45/interpolated/4by4\';
str_o = 'X:\MultiView\YaleB_exteneted\right_45\Results\Our\4by4\phase2\';
str_dummy = 'X:\MultiView\YaleB_exteneted\right_45\Results\Our\phase1\4by4\';  %Initial estiamted results
input1 = 'X:/MultiView/YaleB_exteneted/right_45/interpolated/4by4\';
input2 = 'X:/MultiView/YaleB_exteneted/right_45/phase1\'; %Ground truth
if ~exist(str_o, 'dir')
    mkdir(str_o);
end
ext = '.png';
% str_ref_load= 'J:\TraingSet\CMU_Database\Testing\Frontal_4\';%for calcualte the errors
% fname=flist(kk).name;
for ii = 2:27 % 

    str_input= [str_i sprintf('%02d', ii) ext];    
    %strcat(str_i,num2str(ii),'.bmp');
    str_out = [str_o sprintf('%02d', ii) ext];
    %strcat(str_o,num2str(ii),'.bmp');

    str_in_dummy = [str_dummy sprintf('%02d', ii) ext];
    %strcat(str_dummy,num2str(ii),'.bmp');

    str_ref = [input2 sprintf('%02d', ii) ext];
    %strcat(input2,num2str(ii),'.bmp');
    % img1=imread(str_input);
    img_dummy=imread(str_in_dummy);
    img_dummy=double(img_dummy);
    %%%%%
    img1_1=imread(str_input);  %Input testing face
    img1_1=double(img1_1);
    %%%%%
    % Normalization

    img1=double(img_dummy);

    tic;
    %  Illumination normalizaton using Xie's method 
    %
    [lm,lv]=local_mean_var(img1,25);
    toc;
    img1=img1-lm;
    img1=img1./lv;
    %}


    % Find the reference face examples in the gallery set. The path for the
    % gallery set is set in the function "findneighbor_pca". img2, img2hr are
    % the found interpolated LR and HR reference face examples. m,n are the size of
    % the HR image.
    num_neighbor = 5;
    [img2,img2hr,m,n] = findneighbor_pca(img1_1, num_neighbor,input1,input2,1);


    % Denormalization
    %
    img1=img1.*lv;
    img1=img1+lm;
    img2=img2.*repmat(lv(:),[1,size(img2,2)]);
    img2=img2+repmat(lm(:),[1,size(img2,2)]);
    img2hr=img2hr.*repmat(lv(:),[1,size(img2,2)]);
    img2hr=img2hr+repmat(lm(:),[1,size(img2,2)]);
    %}

    img2=uint8(img2);
    img2hr=uint8(img2hr);

    tmp=[];
    tmphr=[];
    ermap=ones(m,n,size(img2,2));
    for i=1:size(img2,2)
        tmp(:,:,i)=reshape(img2(:,i),m,n);
        % tmphr(:,:,i)=reshape(img2hr(:,i),m*mag,n*mag);
        tmphr(:,:,i)=reshape(img2hr(:,i),m,n);
    end


    % Begin to warp the reference images;
    iii=repmat(img1_1,[1 1 3]);
    img1_1=iii;
    img1_1=uint8(img1_1);
    tic;
    %
    for p=1:size(tmp,3)
        %
        iii=tmp(:,:,p);
        iii=repmat(iii,[1 1 3]);
        iii=uint8(iii);
        iiihr=tmphr(:,:,p);
        iiihr=repmat(iiihr,[1 1 3]);
        iiihr=uint8(iiihr);
        [u,v,test_L,error_map(:,:,p)] = optic_flow_me(img1_1,iii); % Generating the flow field [u, v]

    %     str_out_1=strcat(str_o_1,num2str(p),'.bmp');
    %     str_out_2=strcat(str_o_2,num2str(p),'.bmp');
        test_H = uint8( mywarp_rgb( double(iiihr), u, v ) ) ; % Warp the

        img2_out(:,:,p)=test_L(:,:,1);
        img2hr_out(:,:,p)=test_H(:,:,1);  

        imspace(:,p)=reshape(img2_out(:,:,p),m*n,1);  
        imhrspace(:,p)=reshape(img2hr_out(:,:,p),m*n,1);

    %     iii=repmat(iii,[1 1 3]);



        %  
    end
    imspace=double(imspace);
    imhrspace=double(imhrspace);
    ermap=double(ermap);
    toc;


%%%%%%%%%%%%%%%%%%  PCA RECONSTRUCTION
% psi=mean(imspace,2);
% L=imspace-repmat(psi,[1,size(imspace,2)]);
% 
% psi_1=mean(imhrspace,2);
% L_1=imhrspace-repmat(psi_1,[1,size(imhrspace,2)]);
% 
% % psi_e_n=mean(e_nonf,2);
% % L_e_n=e_nonf; %-repmat(psi_e_n,[1,size(e_nonf,2)]);
% 
% % psi_e_f=mean(e_f,2);
% % L_e_f=e_f; %-repmat(psi_e_f,[1,size(e_f,2)]);
% 
% img1=double(img1);
% img1=reshape(img1(:,:,1),m*n,1);
% img1=img1-psi;
% 
% [V,D]=eig(L'*L);
% diagonal=diag(D);
% [tmp,ind]=sort(diagonal,'descend');
% 
% dim=9;               % The number of eigenvectors to be used
% 
% V=V(:,ind(1:dim));
% D=diag(diagonal(ind(1:dim)));
% diagonal=diag(D);
% diagonal=diagonal.^-0.5;
% D2=diag(diagonal);
% 
% 
% E=L*V*D2;
% imspacew=E'*L;
% testw=E'*img1;
% weight=V*D2*testw;

%%%%%%%%%%%%%% error compensation

% e=abs(double(img1)-L*weight-psi);
% [V_e,D_e]=eig(L_e_n'*L_e_n);
% diagonal=diag(D_e);
% [tmp,ind]=sort(diagonal,'descend');
% 
% dim=6;               % The number of eigenvectors to be used
% 
% V_e=V_e(:,ind(1:dim));
% D_e=diag(diagonal(ind(1:dim)));
% diagonal=diag(D_e);
% diagonal=diagonal.^-0.5;
% D2_e=diag(diagonal);
% 
% E=L_e_n*V_e*D2_e;
% imspacew=E'*L_e_nonf;
% testw=E'*e;
% weight_e=V_e*D2_e*testw;
% e_total=L_e_f*weight_e;

% e_total=reshape(e_total,m*n,1);
% imshow(uint8(e_total));
% e_total=uint8(e_total);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R1=imspacew*imhrspace';%+e_total;
% R1=reshape(R1,m,n);

% R1=double(R1)+double(e_total);

% imshow(uint8(R1))
% imwrite(uint8(R1),str_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refine the structure%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    dx=3; % define the half size of the square neighborhood window
    dy=3;
    sigma=sigma2;
    w=zeros(m,n,2*dx+1,2*dy+1);
    img_dummy=reshape(img_dummy,m*n,1);
    img_dummy=double(img_dummy);
    for i=1:size(imhrspace,2)
        err(:,i)=abs(imhrspace(:,i)-img_dummy);
        tmphr(:,:,i)=reshape(imhrspace(:,i),m,n);
        ermap_t(:,:,i)=reshape(err(:,i),m,n);
    end

%
    tic;
    [B,C]=wcg(tmphr,sigma,2,1,ermap_t);  % learn the local pixel structures;

    w(2:end-1,2:end-1,1,1)=B{1,1};
    w(2:end-1,2:end-1,1,2)=B{1,2};
    w(2:end-1,2:end-1,1,3)=B{1,3};
    w(2:end-1,2:end-1,2,1)=B{2,1};
    w(2:end-1,2:end-1,2,2)=B{2,2};
    w(2:end-1,2:end-1,2,3)=B{2,3};
    w(2:end-1,2:end-1,3,1)=B{3,1};
    w(2:end-1,2:end-1,3,2)=B{3,2};
    w(2:end-1,2:end-1,3,3)=B{3,3};
    toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    img_dummy=reshape(img_dummy,m,n);
    img1hr=double(img_dummy); 
    iter_array= 30;%[1 10 20 30 40 50 60 70 80 90 100];
    for i=1:length(iter_array)
        iter=iter_array(i);  %20  %10
        sd2=0;
        sd=0;

        ww=cell(3,3);
        ww{1,1}=w(:,:,1,1);
        ww{1,1}(1:end-1,1:end-1)=ww{1,1}(2:end,2:end);
        ww{1,2}=w(:,:,1,2);
        ww{1,2}(1:end-1,1:end)=ww{1,2}(2:end,1:end);
        ww{1,3}=w(:,:,1,3);
        ww{1,3}(1:end-1,2:end)=ww{1,3}(2:end,1:end-1);
        ww{2,1}=w(:,:,2,1);
        ww{2,1}(1:end,1:end-1)=ww{2,1}(1:end,2:end);
        ww{2,2}=w(:,:,2,2);
        ww{2,3}=w(:,:,2,3);
        ww{2,3}(1:end,2:end)=ww{2,3}(1:end,1:end-1);
        ww{3,1}=w(:,:,3,1);
        ww{3,1}(2:end,1:end-1)=ww{3,1}(1:end-1,2:end);
        ww{3,2}=w(:,:,3,2);
        ww{3,2}(2:end,1:end)=ww{3,2}(1:end-1,1:end);
        ww{3,3}=w(:,:,3,3);
        ww{3,3}(2:end,2:end)=ww{3,3}(1:end-1,1:end-1);

        tic;
        tol=.01;
        for mm=1:iter
            sd2=sd;
            sd=0;
            ss=zeros(m-2,n-2);
            tt=ww{1,1}.*img1hr;
            ss=ss+tt(1:end-2,1:end-2);
            tt=ww{1,2}.*img1hr;
            ss=ss+tt(1:end-2,2:end-1);
            tt=ww{1,3}.*img1hr;
            ss=ss+tt(1:end-2,3:end);
            tt=ww{2,1}.*img1hr;
            ss=ss+tt(2:end-1,1:end-2);
            tt=ww{2,2}.*img1hr;
            ss=ss+tt(2:end-1,2:end-1);
            tt=ww{2,3}.*img1hr;
            ss=ss+tt(2:end-1,3:end);
            tt=ww{3,1}.*img1hr;
            ss=ss+tt(3:end,1:end-2);
            tt=ww{3,2}.*img1hr;
            ss=ss+tt(3:end,2:end-1);
            tt=ww{3,3}.*img1hr;
            ss=ss+tt(3:end,3:end);
            dd=ss-img1hr(2:end-1,2:end-1);

            img1hr(2:end-1,2:end-1)=img1hr(2:end-1,2:end-1)+0.02*(dd); %-0.3*dd0;    %0.3, 0.1  %%0.055
            sd=sum(sum(abs(dd)));


        % Range constraint
        idi=img1hr>255;
        img1hr(idi)=255;
        idi=img1hr<0;
        img1hr(idi)=0;
        % img1hr(1:mag:end,1:mag:end)=double(ini(1:mag:end,1:mag:end)); % Reset the observation points

        if(abs(sd-sd2)<tol)
            break;
        end

        end %iteration ends
    toc;
%     iter_error_1(i,ii) = calculate_error(uint8(img1hr),'MSE',str_ref);
%     iter_error_2(i,ii) = calculate_error(uint8(img1hr),'MSSIM',str_ref);
    % str_iter=strcat('sub_%d_iteration_',num2str(iter_array(i)),'.png');
    % savename = sprintf(str_iter,ii);
    % 
    % imwrite(uint8(img1hr),fullfile(str_o,savename));
    % clear str_iter
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % imshow(uint8(img1hr))
    if ~exist(str_o)
        mkdir(str_o)
    end
    imwrite(uint8(img1hr),str_out);
    fprintf('N=%d\n',ii);
%     savename_measurement=fullfile(str_o,sprintf('iter_mse+ssim.mat'));
%     save(savename_measurement,'iter_error_1','iter_error_2');
end
% error=calculate_error(img1hr,'MSE','J:\TraingSet\CMU_Database\Testing\Frontal\1.bmp')