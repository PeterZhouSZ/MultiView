function demoMotionVector()
Pts = 3.*rand(3,3);
plot3(Pts(:,1),Pts(:,2),Pts(:,3), 'ro', 'MarkerSize', 20, 'LineWidth', 3); grid on; hold on;
plot3(Pts(1:2,1),Pts(1:2,2),Pts(1:2,3), 'g', 'LineWidth', 3);
plot3(Pts(2:3,1),Pts(2:3,2),Pts(2:3,3), 'g', 'LineWidth', 3);
plot3(Pts([1 3],1),Pts([1 3],2),Pts([1 3],3), 'g', 'LineWidth', 3);
R = eye(3);%[1 0.002 0.0003; -0.0001 1 0.00004; 0 0.00003 1];
T1 = ones(3,1);
T1 = T1 + 1.2*rand(3,1);
newPts = R*Pts + repmat(T1, 1, size(Pts,2));
plot3(newPts(:,1),newPts(:,2),newPts(:,3), 'bo', 'MarkerSize', 20,'LineWidth', 3); grid on; hold on;
plot3(newPts(1:2,1),newPts(1:2,2),newPts(1:2,3), 'g', 'LineWidth', 3);
plot3(newPts(2:3,1),newPts(2:3,2),newPts(2:3,3), 'g', 'LineWidth', 3);
plot3(newPts([1 3],1),newPts([1 3],2),newPts([1 3],3), 'g', 'LineWidth', 3);

T2 = -1*ones(3,1);
T2 = T2 - 1.2*rand(3,1);
newPts1 = R*Pts + repmat(T2, 1, size(Pts,2));
plot3(newPts1(:,1),newPts1(:,2),newPts1(:,3), 'mo', 'MarkerSize', 20, 'LineWidth', 3); grid on; hold on;
plot3(newPts1(1:2,1),newPts1(1:2,2),newPts1(1:2,3), 'g', 'LineWidth', 3);
plot3(newPts1(2:3,1),newPts1(2:3,2),newPts1(2:3,3), 'g', 'LineWidth', 3);
plot3(newPts1([1 3],1),newPts1([1 3],2),newPts1([1 3],3), 'g', 'LineWidth', 3);

end

% function drawPts(Pts)
%     plot3(Pts(:,1),Pts(:,2),Pts(:,3), 'ro', 'MarkerSize', 20); grid on; hold on;
%     plot3(Pts(1:2,1),Pts(1:2,2),Pts(1:2,3), 'g', 'LineWidth', 3);
%     plot3(Pts(2:3,1),Pts(2:3,2),Pts(2:3,3), 'g', 'LineWidth', 3);
%     plot3(Pts([1 3],1),Pts([1 3],2),Pts([1 3],3), 'g', 'LineWidth', 3);
% end