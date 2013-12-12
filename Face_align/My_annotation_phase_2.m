function My_annotation_phase_2(path)
%
% function annotate( imagetype )
%
%
% Reads and displays all images in the current dir one by one.
%
% When an image is displayed, annotation can be done by left-
% clicking the mouse. Press 'e' to end annotation and write
% a corresponding file in the current dir.
% 
%
% $Id: annotate.m,v 1.1.1.1 2003/01/03 19:18:51 aam Exp $ 
%

imagetype='*.bmp';  % the format of the images in current dir.
%path='E:\HZ\neighbor embedding\Face_align\database\phase_1_adjustment\';

% eval( sprintf('dirlist=dir(''*.%s'');', imagetype ) ); % get a dir list of all images of type 'imagetype'
dirlist=dir([path,imagetype]);
% x=[];
x=struct('name','nnnnnxxfffq_yymmdd','Leyex',-1,'Leyey',-1,'Reyex',-1,'Reyey',-1,'Cmouthx',-1,'Cmouthy',-1);%,'Cnosex',-1,'Cnosey',-1
list=[];
N=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each image in the current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 40%1 :length(dirlist)

	%%%%%%%%%%%%%%%%%%%%%%
	% load and show image
	%%%%%%%%%%%%%%%%%%%%%%      
    image_filename=dirlist(i).name;   
    % if(image_filename([6 7 9])=='fa1')
    % eval(sprintf('[img,cmap]=imread(''%s'');',image_filename));
    [img,cmap]=imread([path,image_filename]);
    % figure;    
	axis ij;	
    imagesc(img);
    % imshow(img,'InitialMagnification',200)
    axis image;
    channels = size(img,3);
    if channels==1,    
        if isempty(cmap),       
            colormap( gray );
        else
            colormap( cmap );
        end         
    end
    eval(sprintf('title(''%s'');',image_filename));
    s=size(img);
    height=s(1);
    width=s(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % annotate image with points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end_botton = 0;
    point_count = 0;
    path_count = 0;
    while end_botton~='e' && end_botton~='b'
    
    	[points, end_botton]=markpath( point_count );
    	s=size(points);            	    
    	npathpoints = s(1);    	
    	if ((npathpoints==0)||(npathpoints>3))    %%%%%%%%%%%%Take Care 	    	
    		break;
		end    		
    end
    
    if(((npathpoints==0)||(npathpoints>3)) ) %%%%Take care of numbers!!!
        continue;
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Construct the list.mat
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    x.name=image_filename;
    
    x.Leyex=uint16(points(1,1));
    x.Leyey=uint16(points(1,2));
    x.Reyex=uint16(points(2,1));
    x.Reyey=uint16(points(2,2));

    
    x.Cmouthx=uint16(points(3,1));
    x.Cmouthy=uint16(points(3,2));
    
%     x.Cnosex=uint16(points(4,1));
%     x.Cnosey=uint16(points(4,2));
%     x(1,1)=uint16(points(1,1));
%     x(2,1)=uint16(points(1,2));
%     x(1,2)=uint16(points(2,1));
%     x(2,2)=uint16(points(2,2));
    list=[list;x];
%     save([path,image_filename(1:end-4)],'x')
    N=N+1;
    fprintf('N=%d\n',N);

    end
   
    
    
end
save([path,'Myannotate.mat'],'list')



