clc
clear all
close all

set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',12)

%% load data
saveFigFlagOuter = 0;
labelsize = 14;
legendsize = 14;

load('./Data/ROAData.mat')



% Plot positive points
figure(1)
figSize = [0 0 0.35 0.2]*.75;
set(gcf,'units','normalized','position',figSize)
count_pos = length(XMatrix_pos(1,:));
count_neg = length(XMatrix_neg(1,:));
for i = 1:count_pos
    X_i = XMatrix_pos(:,i);
    plot(X_i(1),X_i(2),'g.','MarkerSize',1);
    if i == 1
        hold on; box on; grid on;
    end
end

% Plot negative points
for i = 1:count_neg
    X_i = XMatrix_neg(:,i);
    plot(X_i(1),X_i(2),'r.','MarkerSize',1);
end
% plot(specialX(1),specialX(2),'k.','markersize',15)

% Plot ghost points for the legend
h2 = plot(nan,nan,'g.','MarkerSize',10);
h1 = plot(nan,nan,'r.','MarkerSize',10);

% Plot sublevel set
VBar = 0.8*min(costNegVec);
[~,h3] = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-');
legend([h3 h2 h1],'$\Upsilon_N$','$x \in \Gamma_N$','$x \notin \Gamma_N$','interpreter','Latex','Fontsize',legendsize,'location','northeast')
xlabel('$\tilde{x}_1$','interpreter','Latex','Fontsize',labelsize)
ylabel('$\tilde{x}_2$','interpreter','Latex','FontSize',labelsize)

if saveFigFlagOuter == 1
    figure(1)
    filename = strcat('./Plots/ROAPlot_doubleInt');
    saveas(gcf,filename,'epsc')
end
%}

%% Plot phase plane
colorMatrix = colororder;
greyColor = [111 111 111]/255;

figSize = [0 0 0.35 0.5];
subPlotGap = [0.08 0.1];
subPlotH = [0.1 0.04];
subPlotW = [0.1 0.05];

markerStr = cell(4,1);
markerStr{1} = 'o';
markerStr{2} = 'd';
markerStr{3} = '^';
markerStr{4} = 's';
sillyMarkerSize = 6;

% Plot the sublevel set
VBar = 0.8*min(costNegVec);
figure(2)
set(gcf,'units','normalized','position',figSize)
subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
[~,h0] = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-');
hold on; grid on; box on;
subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW)
[~,h0_2] = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',[0 0 0],'linestyle','-');
hold on; grid on; box on;
for iOuter = 1:4
    % Load tracking
    loadstr = ['./Data/output_tracking_N10_ic',num2str(iOuter)];
    load(loadstr)
    output_tracking = output;

    % Load IG
    loadstr = ['./Data/output_ig_N10_ic',num2str(iOuter),'_ws0'];
    load(loadstr)
    output_ig = output;

    % Plot trajectories
    Gx = output_tracking.controlArgs.Gx;
    r = output_tracking.controlArgs.r(1);
    X_tracking = output_tracking.X*0;
    X_ig = output_ig.X*0;
    output_ig.refHist(end) = output_ig.refHist(end-1);
    for i = 1:length(output_tracking.X(1,:))
        X_tracking(:,i) = output_tracking.X(:,i) - Gx*r;
        X_ig(:,i) = output_ig.X(:,i) - Gx*output_ig.refHist(i);
    end
    figure(2)
    subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
    hold on;
    h1 = plot(X_tracking(1,:),X_tracking(2,:),'linewidth',1,'Color',colorMatrix(1,:));
    h4 = plot(X_tracking(1,1),X_tracking(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(1,:),'linewidth',1);
    plot(X_tracking(1,1),X_tracking(2,1),'.','MarkerSize',5,'Color',colorMatrix(1,:),'linewidth',1)
    plot(0,0,'k.','MarkerSize',10);

    % % Plot also in x plane
    h2 = plot(output_ig.X(1,:),output_ig.X(2,:),'Color',colorMatrix(4,:),'linewidth',1);

    %     h3 = plot(X_ig(1,:),X_ig(2,:),'-','linewidth',1,'Color',colorMatrix(2,:));
    % h5 = plot(X_ig(1,1),X_ig(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(2,:),'linewidth',1);
    % plot(X_ig(1,1),X_ig(2,1),'.','MarkerSize',5,'Color',colorMatrix(2,:),'linewidth',1)


    subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW)
    h3 = plot(X_ig(1,:),X_ig(2,:),'-','linewidth',1,'Color',colorMatrix(2,:));
    plot(X_ig(1,1),X_ig(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(2,:),'linewidth',1);
    plot(X_ig(1,1),X_ig(2,1),'.','MarkerSize',5,'Color',colorMatrix(2,:),'linewidth',1)
end

subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
xlim([-1.1 1.1])
ylim([-0.25 0.25])
legend([h0 h1 h2],'$\psi(x) = r_\psi$','$x_k$ (rMPC)','$x_k$ (IG-MPC)','interpreter','Latex','Fontsize',legendsize,'numcolumns',2)
xlabel('$x_1$','interpreter','Latex','Fontsize',labelsize)
ylabel('$x_2$','interpreter','Latex','FontSize',labelsize)


subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW)
xlim([-0.12 0.12])
ylim([-0.1 0.1])
plot(0,0,'k.','MarkerSize',10);
xlabel('$\tilde{x}_1$','interpreter','Latex','Fontsize',labelsize)
ylabel('$\tilde{x}_2$','interpreter','Latex','FontSize',labelsize)
legend([h0_2 h3],'$\psi(\tilde{x}) = r_\psi$','$\tilde{x}_k$ (IG-MPC)','interpreter','Latex','Fontsize',legendsize,'numcolumns',2)
% legend(h3,'$\tilde{x}_k$ (IG-MPC)','interpreter','Latex','Fontsize',legendsize,'location','northeast', 'numcolumns',2)

if saveFigFlagOuter == 1
    figure(2)
    filename = strcat('./Plots/PhasePlot_doubleInt');
    saveas(gcf,filename,'epsc')
end
%}


%{
%% Plot trajectories
figSize = [0 0 0.35 0.25];
subPlotGap = [0.1 0.1];
subPlotH = [0.22 0.04];
subPlotW = [0.1 0.05];
sillyMarkerSize = 6;


figure(3)
    set(gcf,'units','normalized','position',figSize)
% s
for iOuter = 1:4
    % Load IG
    loadstr = ['output_ig_N10_ic',num2str(iOuter),'_ws1'];
    load(loadstr)

    t = output.t;
    umax = output.controlArgs.umax;
    umin = output.controlArgs.umin;
    X = output.X;
    U = output.U;
    refHist = output.refHist;
    ellHist = output.ellHist;

    % x1
    subtightplot(2,3,1,subPlotGap,subPlotH,subPlotW)
    plot(t,X(1,:),'Color',colorMatrix(iOuter,:))
    hold on; box on; grid on
    plot(t,refHist,'Linestyle','-.','LineWidth',1,'Color',colorMatrix(iOuter,:))
    plot(t(1),X(1,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    plot(t(1),X(1,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % x2
    subtightplot(2,3,2,subPlotGap,subPlotH,subPlotW)
    plot(t,X(2,:),'Color',colorMatrix(iOuter,:))
        hold on; box on; grid on
    plot(t(1),X(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    plot(t(1),X(2,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % delta
    subtightplot(2,3,3,subPlotGap,subPlotH,subPlotW)
    if iOuter == 1
        h_u = plot(t,t*0+umax(1),'k-.','LineWidth',1);
        hold on; box on; grid on;
        plot(t,t*0+umin(1),'k-.','LineWidth',1)
    end
    plot(t(1:end-1),U(1,:),'LineWidth',1,'Color',colorMatrix(iOuter,:))

     % VHist
    subtightplot(2,3,4,subPlotGap,subPlotH,subPlotW);
    if iOuter == 1
        h_psi = plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.r_psi,'k-.');
        hold on; box on; grid on;
    end
    plot(t(1:end-1),sqrt(output.VHist),'LineWidth',1,'color',colorMatrix(iOuter,:));

     % Error
    subtightplot(2,3,5,subPlotGap,subPlotH,subPlotW);
    if iOuter == 1
        h_phi = semilogy(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.r_phi,'k-.');
        hold on; box on; grid on;
    end
        semilogy(t(1:end-1),output.errorHist,'LineWidth',1,'color',colorMatrix(iOuter,:));

     % Iterations
    subtightplot(2,3,6,subPlotGap,subPlotH,subPlotW);
    if iOuter == 1
         h_ell = plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.ellStar_ws,'k-.');
         hold on; box on; grid on;
    end
    stairs(t(1:end-1),ellHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
end

% Add legends and axes
% x1
subtightplot(2,3,1,subPlotGap,subPlotH,subPlotW)
ylabel('$x_1$','interpreter','Latex','FontSize',labelsize)
% legend('$x_{1,k}$','$v_k$','interpreter','Latex','Fontsize',legendsize,'location','northeast')
ylim([-1.1 1.1])

% x2
subtightplot(2,3,2,subPlotGap,subPlotH,subPlotW)
ylabel('$x_2$','interpreter','Latex','FontSize',labelsize)
ylim([-0.1 0.1])

% u 
subtightplot(2,3,3,subPlotGap,subPlotH,subPlotW)
% legend(h_u,'$\mathcal{U}$','interpreter','Latex','Fontsize',legendsize,'location','southeast')
% xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$u$','interpreter','Latex','FontSize',labelsize)
ylim([-0.11 0.11])

% VHist
subtightplot(2,3,4,subPlotGap,subPlotH,subPlotW);
legend(h_psi,'$r_\psi$','interpreter','Latex','Fontsize',legendsize,'location','southwest')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\psi(\tilde{x})$','interpreter','Latex','FontSize',labelsize)
ylim([0 1.05*output.controlArgs.const.r_psi])

% Error hist
subtightplot(2,3,5,subPlotGap,subPlotH,subPlotW);
legend(h_phi,'$r_\phi$','interpreter','Latex','Fontsize',legendsize,'location','southwest')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\phi(\tilde{x},\hat{\textnormal{u}})$','interpreter','Latex','FontSize',labelsize)
ylim([1e-12 1e-2])
yticks([1e-12 1e-9 1e-6 1e-3])

% Iterations
subtightplot(2,3,6,subPlotGap,subPlotH,subPlotW);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell$','interpreter','Latex','FontSize',labelsize)
ylim([0 35])
yticks([0 10 20 30])
legend(h_ell,'$\ell^*_{\textrm{ws}}$','interpreter','Latex','Fontsize',legendsize,'location','northeast')


if saveFigFlagOuter == 1
    figure(3)
    filename = strcat('./Plots/wsTraj_doubleInt');
    saveas(gcf,filename,'epsc')
end

%}

%% Plot trajectories (CS)
figSize = [0 0 0.35 0.5];
subPlotGap = [0.08 0.1];
subPlotH = [0.1 0.04];
subPlotW = [0.1 0.05];
sillyMarkerSize = 6;

figure(4)
    set(gcf,'units','normalized','position',figSize)
% s
for iOuter = 1:4
    % Load IG
    loadstr = ['./Data/output_ig_N10_ic',num2str(iOuter),'_ws0'];
    load(loadstr)

   
    t = output.t;
    umax = output.controlArgs.umax;
    umin = output.controlArgs.umin;
    X = output.X;
    U = output.U;
    refHist = output.refHist;
    ellHist = output.ellHist;

    % x1
    subtightplot(3,2,1,subPlotGap,subPlotH,subPlotW)
    plot(t,X(1,:),'Color',colorMatrix(iOuter,:))
    hold on; box on; grid on
    plot(t,refHist,'Linestyle','-.','LineWidth',1,'Color',colorMatrix(iOuter,:))
    plot(t(1),X(1,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    plot(t(1),X(1,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % x2
    subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
    plot(t,X(2,:),'Color',colorMatrix(iOuter,:))
        hold on; box on; grid on
    plot(t(1),X(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    plot(t(1),X(2,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % delta
    subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
    if iOuter == 1
        h_u = plot(t,t*0+umax(1),'k-.','LineWidth',1);
        hold on; box on; grid on;
        plot(t,t*0+umin(1),'k-.','LineWidth',1)
    end
    plot(t(1:end-1),U(1,:),'LineWidth',1,'Color',colorMatrix(iOuter,:))

     % VHist
    subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
    plot(t(1:end-1),sqrt(output.VHist)/output.controlArgs.const.r_psi,'LineWidth',1,'color',colorMatrix(iOuter,:));
    hold on; box on; grid on;

     % Error
    subtightplot(3,2,5,subPlotGap,subPlotH,subPlotW);
    semilogy(t(1:end-1),output.errorHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
    hold on; box on; grid on;

     % Iterations
    subtightplot(3,2,6,subPlotGap,subPlotH,subPlotW);
    if iOuter == 1
         h_ell = plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.ellStar_cs,'k-.');
         hold on; box on; grid on;
    end
    stairs(t(1:end-1),ellHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
end
% Add legends and axes
% x1
subtightplot(3,2,1,subPlotGap,subPlotH,subPlotW)
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$x_1$','interpreter','Latex','FontSize',labelsize)
% legend('$x_{1,k}$','$v_k$','interpreter','Latex','Fontsize',legendsize,'location','northeast')
ylim([-1.1 1.1])

% x2
subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$x_2$','interpreter','Latex','FontSize',labelsize)
ylim([-0.1 0.1])

% u 
subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
% legend(h_u,'$\mathcal{U}$','interpreter','Latex','Fontsize',legendsize,'location','southeast')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$u$','interpreter','Latex','FontSize',labelsize)
ylim([-0.105 0.105])

% VHist
subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
% legend(h_psi,'$r_\psi$','interpreter','Latex','Fontsize',legendsize,'location','southwest')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\psi(\tilde{x})/r_\psi$','interpreter','Latex','FontSize',labelsize)
ylim([0 1.05])

% Error hist
subtightplot(3,2,5,subPlotGap,subPlotH,subPlotW);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\phi(\tilde{x},\hat{\textnormal{u}})$','interpreter','Latex','FontSize',labelsize)
ylim([-Inf 1e-1])
yticks([1e-5 1e-3 1e-1])

% Iterations
subtightplot(3,2,6,subPlotGap,subPlotH,subPlotW);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell$','interpreter','Latex','FontSize',labelsize)
ylim([0 35])
yticks([0 10 20 30])
legend(h_ell,'$\ell^*$','interpreter','Latex','Fontsize',legendsize,'location','northeast')

if saveFigFlagOuter == 1
    figure(4)
    filename = strcat('./Plots/csTraj_doubleInt');
    saveas(gcf,filename,'epsc')
end

%% Plot the nested ROAs
NVec = [10 25 50 75];
figure(5)

figSize = [0 0 0.35 0.5];
subPlotGap = [0.08 0.1];
subPlotH = [0.1 0.04];
subPlotW = [0.1 0.05];

set(gcf,'units','normalized','position',figSize)
legVec = cell(length(NVec),1);
for iOuter = 1:length(NVec)
    loadstr = ['./Data/ROAData_N',num2str(NVec(iOuter))];
    load(loadstr);
    subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
    VBar = 0.8*min(costNegVec);
    [~,h0] = contour(XGrid,YGrid,costGrid',[VBar VBar],'linewidth',1,'color',colorMatrix(iOuter,:),'linestyle','-');
    hold on; grid on; box on;
    legVec{iOuter} = ['$N = ',num2str(NVec(iOuter)),'$'];
end
xlabel('${x}_1$','interpreter','Latex','Fontsize',labelsize)
ylabel('${x}_2$','interpreter','Latex','FontSize',labelsize)
legend(legVec,'interpreter','Latex','Fontsize',legendsize,'location','northeast','orientation','horizontal')

% Plot ellStar trend
load('./Data/ellStarTrend.mat')
% subPlotH = [0.22 0.04];
subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW)
iPlot = 40;
% semilogy(NVec(1:iPlot),ellStar_ws(1:iPlot),'color',colorMatrix(1,:))
semilogy(NVec(1:iPlot),ellStar_cs(1:iPlot),'k')
hold on; box on; grid on;
%  legend('$\ell^*_{\textrm{ws}}$','$\ell^*_{\textrm{cs}}$','interpreter','Latex','Fontsize',legendsize,'location','northwest')
% xlim([2 Inf])
xlabel('$N$','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell^*$','interpreter','Latex','FontSize',labelsize)
% semilogy(NVec(1:iPlot),ellStar_cs(1:iPlot)-ellStar_ws(1:iPlot),'color',[0 0 0])
xlim([5 40])
yticks([1 1e1 1e2 1e3 1e4 1e5])
% 
% % PLOT MINI PLOT
% scale = 0.6;
% ax2 = axes('Position',[.75 .3 0.3*scale 0.5*scale]);
% iPlot = 40;
% semilogy(NVec(1:iPlot),ellStar_cs(1:iPlot)-ellStar_ws(1:iPlot),'color',[0 0 0])
% hold on; box on; grid on;
% xlim([2 Inf])
% % ylimmini = [-2 4];
% % xlimmini = [20 30];
% % xlim(xlimmini)
% % ylim(ylimmini)
% set(gca,'YTick',[1 10 100],'XTick',[2 40])
% % set(gca,'XTick',[2 4 6 8 10])
% % ylim([1 Inf])
% ylabel('$\ell^*_{\textrm{cs}} - \ell^*_{\textrm{ws}}$','interpreter','Latex','FontSize',labelsize)

if saveFigFlagOuter == 1
    figure(5)
    filename = strcat('./Plots/NTrend_doubleInt');
    saveas(gcf,filename,'epsc')
end
%}