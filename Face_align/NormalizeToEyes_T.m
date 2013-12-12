function A = NormalizeToEyes_T(img, LeyeX, LeyeY, ReyeX, ReyeY, r1, r2, r3, newW, newH)
%NormalizeToEyes_T      Normalize face image to the eyes
% Refer to the paper, K.-W. Wong, K.-M. Lam and W.-C. Siu, An efficient algorithm for human face detection and facial feature extraction under different conditions, Pattern recognition 34, pp. 1993-2004, 2001.
% A = NormalizeToEyes_T(img, LeyeX, LeyeY, ReyeX, ReyeY, r1, r2, r3, newW, newH)
%   A = output image in 2D matrix
%   img = input image in 2D matrix
%   LeyeX = x-coordinate of eye on left side of image
%   LeyeY = y-coordinate of eye on left side of image
%   ReyeX = x-coordinate of eye on right side of image
%   ReyeY = y-coordinate of eye on right side of image
%   r1 = first ratio, 1.8, equation 4a
%   r2 = second ratio, 0.2, equation 4b
%   r3 = third ratio, 0.3, equation 4c
%   newW = new width of output image
%   newH = new height of output image
% by Thomas Tse
% Last updated on 27 March 2007

% on input image
[imgH,imgW] = size(img);
a = ReyeY - LeyeY;
b = ReyeX - LeyeX;
dEye = sqrt(a^2 + b^2);
phi = atan2(a,b);
% on output image
hFace = round(r1 * dEye);
hEye = round(r2 * hFace);
wEye = round(r3 * hFace);
wFace = round(wEye + wEye + dEye);

% for fast computing in looping
hEye_p1 = hEye + 1;
wEye_p1 = wEye + 1;
% cosp = cos(phi);
% sinp = sin(phi);
% LeyeX_p1 = LeyeX + 1;
% LeyeY_p1 = LeyeY + 1;
Mat =[cos(phi), -sin(phi); sin(phi), cos(phi)];
Leye_p1 = [LeyeX;LeyeY]+1;
x = 1:wFace;
xx = x - wEye - 1;% shift the centre to eye

% B = uint8(zeros(hFace,wFace));% initialize
for y = 1:hFace
%     % faster, but has error, don't know why
%     y1 = y .* ones(1,wFace);
%     yy = y1 - hEye_p1;
%     x = 1:wFace;
%     % shift the centre to eye
%     xx = x - wEye_p1;
%     % rotate at the eye
%     idxXr = round( xx*cosp - yy*sinp);
%     idxYr = round( xx*sinp + yy*cosp);
%     idxXr = idxXr + LeyeX_p1;
%     idxYr = idxYr + LeyeY_p1;
%     % check
%     ii = find(idxXr > 1 & idxXr<= imgW & idxYr > 1 & idxYr <= imgH);
%     B(y1(ii),x(ii)) = img(idxYr(ii),idxXr(ii));

    % faster, have just modified, no error
    yy = y .* ones(1,wFace) - hEye_p1;% shift the centre to eye
    % rotate at the eye
    idxr = Mat*[xx;yy];% rotate the index plane
    idxr = idxr + Leye_p1*ones(1,wFace);% adjust to the left eye coordinates
    idxr = round(idxr);% make it integer
    % check
    ii = find(idxr(1,:) > 1 & idxr(1,:)<= imgW & idxr(2,:) > 1 & idxr(2,:) <= imgH);
    idxrr = [];% initialize the vector
    idxrr(ii) = idxr(2,ii) + (idxr(1,ii)-1) * imgH;% cast to 1D vector (must)
    B(y,x(ii)) = img(idxrr(ii));% copy
    
%     % no error, but slower
%     % put the window at the eye in the image
%     yy = y - hEye_p1;
%     for x = 1:wFace
%         % shift the centre to eye
%         xx = x - wEye_p1;
%         % rotate at the eye
%         idxXr = round( cosp*xx - sinp*yy );
%         idxYr = round( sinp*xx + cosp*yy );
%         % put the window at the eye in the image
%         idxXr = idxXr + LeyeX_p1;
%         idxYr = idxYr + LeyeY_p1;
%         % check if exceed the image
%         if(idxXr > 1 & idxXr<= imgW & idxYr > 1 & idxYr <= imgH),
%             B(y,x) = img(idxYr,idxXr);
%         else, B(y,x) = 0;
%         end
%     end
    
%     % no error, but slower
%     % put the window at the eye in the image
%     yy = y - hEye_p1;
%     idxr2=[];
%     for x = 1:wFace
%         % shift the centre to eye
%         xx = x - wEye_p1;
%         % rotate at the eye
%         idxr = floor( Mat*[xx;yy] );
%         % put the window at the eye in the image
%         idxr = idxr + Leye_p1;
% 
%         idxr2=[idxr2,idxr];
%         % check if exceed the image
%         if(idxr(1) > 1 & idxr(1)<= imgW & idxr(2) > 1 & idxr(2) <= imgH),
%             B(y,x) = img(idxr(2),idxr(1));
%         else, B(y,x) = 0;
%         end
%     end

end

A=imresize(B,[newH,newW],'bicubic');
% A = ImageResampling_T(B, newW, newH);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Out = ImageResampling_T(in, new_width, new_height)
% image resampling

[old_height, old_width] = size(in);
ratioW = old_width / new_width;
ratioH = old_height / new_height;
if(ratioW>1),preciW=0.5;else,preciW=0;end
if(ratioH>1),preciH=0.5;else,preciH=0;end

w = 1:new_width;
w1 = floor((w-1)*ratioW+preciW)+1;
if(w1>old_width) 
    w1=old_width;
end
for h = 1:new_height
    h1 = floor((h-1)*ratioH+preciH)+1;
    if(h1>old_height)
        h1=old_height;
    end
    Out(h,w) = in(h1,w1);
end
