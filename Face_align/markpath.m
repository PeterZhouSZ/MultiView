function [points, end_botton]=markpath( point_count )
%
% function [points, end_botton]=markpath( point_count )
%
%
% Annotates one open or closed path and returns the point 
% coordinates and the keyboard press that ended the 
% annotation.
%
% 'point_count' is the current point count before invoking markpath()
%
% Output format: points = [ x y ; x y ; .... ]
%
% $Id: markpath.m,v 1.1.1.1 2003/01/03 19:18:51 aam Exp $ 
%

% setup
v=[];
left_mouse   = 1;
middle_mouse = 2;
right_mouse  = 3;

if point_count==0

	% output help to user
	disp('left-click  : Mark a point');
	disp('right-click : Undo last point');	
	disp('   ''o''      : End this (open) path');
	disp('   ''c''      : End and close this path');	
	disp('   ''e''      : End this shape');	
	disp('');
	disp('NOTICE: Closed paths should be defined clock-wise.');
end

% get userinput
while 1,	   
    [x,y,input]=ginput(1);        	    
    if input==left_mouse 			% add point	   
       v=[v; x y];
       line( v(:,1), v(:,2), 'Marker', '+','Color','red', 'LineWidth', 2 );
    elseif input==right_mouse		% undo point
       line( v(:,1), v(:,2), 'Marker', '+','Color','green', 'LineWidth', 2 );
       v=v(1:end-1,:); 				% remove the last point
       line( v(:,1), v(:,2), 'Marker', '+','Color','red', 'LineWidth', 2 );
    elseif input=='c' 				% close path           
		line( [ v(:,1); v(1,1)], [v(:,2); v(1,2)], 'Marker', '+','Color','red', 'LineWidth', 2 );           
		break;
    elseif (input=='e' | input=='o') % stop annotating this shape        
       break;
    end    
    eval(sprintf('xlabel(''# points=%i'');', point_count+size(v,1)));    
end		
points=v;
end_botton = input;