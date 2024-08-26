clc
clear all
close all

%% load data
load('ROAData.mat')

% Plot positive points
% count_pos = length(XMatrix_pos(1,:));
% count_neg = length(XMatrix_neg(1,:));
% for i = 1:count_pos
%     X_i = XMatrix_pos(:,i);
%     plot(X_i(1),X_i(2),'g.','MarkerSize',1);
%     if i == 1
%         hold on; box on; grid on;
%     end
% end

% % Plot negative points
% for i = 1:count_neg
%     X_i = XMatrix_neg(:,i);
%     plot(X_i(1),X_i(2),'r.','MarkerSize',1);
% end
% plot(specialX(1),specialX(2),'k.','markersize',15)


% Load the batch data and plot
load('batchData_tracking_N10.mat')
figure(1)
hold on
Gx = batchData.controlArgs.Gx;
r = [0; 0];
NSample_batch = length(batchData.output);
for iOuter = 1:NSample_batch
    X_traj = batchData.output(iOuter).X; % XTilde and x are equivalent in this case
    if norm(X_traj(:,end)) > 0.01
        plot(X_traj(1,:),X_traj(2,:),'r-');
        % plot(X_traj(1,1),X_traj(2,1),'r.','MarkerSize',1)
    else
        plot(X_traj(1,:),X_traj(2,:),'g-');
        % plot(X_traj(1,1),X_traj(2,1),'g.','MarkerSize',1)
    end
end
xlim([-1 1])
ylim([-0.3 0.3])

% Plot the sublevel set
VBar = 0.9*min(costNegVec);
h0 = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-.');



% Load the batch data and plot
load('batchData_fg_N10.mat')
figure(2)
hold on
Gx = batchData.controlArgs.Gx;
r = [0; 0];
NSample_batch = length(batchData.output);
for iOuter = 1:NSample_batch
    X_traj = batchData.output(iOuter).X; % XTilde and x are equivalent in this case
    if norm(X_traj(:,end)) > 0.01
        plot(X_traj(1,1),X_traj(2,1),'r.');
        % plot(X_traj(1,1),X_traj(2,1),'r.','MarkerSize',1)
    else
        plot(X_traj(1,1),X_traj(2,1),'g.');
        % plot(X_traj(1,1),X_traj(2,1),'g.','MarkerSize',1)
    end
end
xlim([-1 1])
ylim([-0.3 0.3])

% Plot the sublevel set
VBar = 0.9*min(costNegVec);
h0 = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-.');

figure(3)
hold on
Gx = batchData.controlArgs.Gx;
r = [0; 0];
NSample_batch = length(batchData.output);
for iOuter = 1:NSample_batch
    X_traj = batchData.output(iOuter).X; % XTilde and x are equivalent in this case
    if norm(X_traj(:,end)) > 0.01
        plot(X_traj(1,:),X_traj(2,:),'r');
        % plot(X_traj(1,1),X_traj(2,1),'r.','MarkerSize',1)
    else
        plot(X_traj(1,:),X_traj(2,:),'g');
        % plot(X_traj(1,1),X_traj(2,1),'g.','MarkerSize',1)
    end
end
xlim([-1 1])
ylim([-0.3 0.3])

% Plot the sublevel set
VBar = 0.9*min(costNegVec);
h0 = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-.');

figure(4)
hold on
Gx = batchData.controlArgs.Gx;
r = [0; 0];
NSample_batch = length(batchData.output);
for iOuter = 1:NSample_batch
    X_traj_non = batchData.output(iOuter).X; % XTilde and x are equivalent in this case
    X_traj = X_traj_non*0;
    for j = 1:length(X_traj_non(1,:))
       X_traj(:,j) = X_traj_non(:,j) - Gx*batchData.output(iOuter).refHist(j);
   end
    if norm(X_traj(:,end)) > 0.01
        plot(X_traj(1,:),X_traj(2,:),'r');
        % plot(X_traj(1,1),X_traj(2,1),'r.','MarkerSize',1)
    else
        plot(X_traj(1,:),X_traj(2,:),'g');
        % plot(X_traj(1,1),X_traj(2,1),'g.','MarkerSize',1)
    end
end
xlim([-1 1])
ylim([-0.3 0.3])

% Plot the sublevel set
VBar = 0.9*min(costNegVec);
h0 = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-.');
