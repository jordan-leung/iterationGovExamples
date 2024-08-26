clc
clear all
close all

set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',12)

%% load data
saveFigFlagOuter = 0;
labelsize = 14;
legendsize = 14;
colorMatrix = colororder;
temp = colorMatrix(3,:);
colorMatrix(3,:) = colorMatrix(4,:);
colorMatrix(4,:) = temp;
greyColor = [111 111 111]/255;

% Load case 
ws_flag = 0;
% ws_flag = 1;


%% Plot trajectories

markerStr = cell(4,1);
markerStr{1} = 'o';
markerStr{2} = 'd';
markerStr{3} = '^';
markerStr{4} = 's';
sillyMarkerSize = 6;
NVec = [15 25 35];


figSize = [0 0 0.35 0.5];
subPlotGap = [0.08 0.1];
subPlotH = [0.1 0.04];
subPlotW = [0.1 0.05];


figure(1)
set(gcf,'units','normalized','position',figSize)
% s
for iOuter = 1:3
    % Load IG
    loadstr = ['./Data/output_ig_N',num2str(NVec(iOuter)),'_ic1_ws',num2str(ws_flag)];
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
    hVec(iOuter) = plot(t,X(1,:),'Color',colorMatrix(iOuter,:));
    hold on; box on; grid on
    plot(t,refHist,'Linestyle','-.','LineWidth',1,'Color',colorMatrix(iOuter,:))
    % plot(t(1),X(1,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    % plot(t(1),X(1,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % x2
    subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
    plot(t,X(3,:)*180/pi,'Color',colorMatrix(iOuter,:))
        hold on; box on; grid on
    % plot(t(1),X(2,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
    % plot(t(1),X(2,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % delta
    subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
    % if iOuter == 1
    %     h_u = plot(t,t*0+umax(1),'k-.','LineWidth',1);
    %     hold on; box on; grid on;
    %     plot(t,t*0+umin(1),'k-.','LineWidth',1)
    % end
    hold on; box on; grid on;
    plot(t(1:end-1),U(1,:),'LineWidth',1,'Color',colorMatrix(iOuter,:))

     % VHist
    subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
    plot(t(1:end-1),sqrt(output.VHist)/output.controlArgs.const.r_psi,'LineWidth',1,'color',colorMatrix(iOuter,:));
    hold on; box on; grid on;

     % Error
    subtightplot(3,2,5,subPlotGap,subPlotH,subPlotW);
    if  ws_flag 
        h_phi = semilogy(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.r_phi,'linestyle','-.','color',colorMatrix(iOuter,:));
        hold on; box on; grid on;
    end
        semilogy(t(1:end-1),output.errorHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
        hold on; box on; grid on;

     % Iterations
     subtightplot(3,2,6,subPlotGap,subPlotH,subPlotW);
     semilogy(t(1:end-1),ellHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
     hold on; box on; grid on;
    % plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.ellStar_ws,'color',colorMatrix(iOuter,:),'linestyle','-.');
end

% % Process tracking plots and add
% load('output_tracking_N30_ic1.mat')
% X = output.X;
% U = output.U;
% plotInd = length(t);


% Add legends and axes
% x1
subtightplot(3,2,1,subPlotGap,subPlotH,subPlotW)
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$s$ [m]','interpreter','Latex','FontSize',labelsize)
if ~ws_flag
    legend(hVec,['$N = ',num2str(NVec(1)),'$'],['$N = ',num2str(NVec(2)),'$'],['$N = ',num2str(NVec(3)),'$'],...
        'interpreter','Latex','Fontsize',legendsize,'location','northeast')
end
xlim([0,10])
% plot(t(1:plotInd),X(1,1:plotInd),'k-')
% ylim([-0.1 1.1])

% x2
subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$\theta$ [deg]','interpreter','Latex','FontSize',labelsize)
xlim([0,10])
% plot(t(1:plotInd),X(3,1:plotInd)*180/pi,'k-')

% u 
subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW)
% legend(h_u,'$\mathcal{U}$','interpreter','Latex','Fontsize',legendsize,'location','southeast')
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$u$ [N]','interpreter','Latex','FontSize',labelsize)
plot(t,t*0+umax(1),'k-.','LineWidth',1);
plot(t,t*0+umin(1),'k-.','LineWidth',1);
ylim([1.05*umin 1.05*umax])
xlim([0,10])


% VHist
subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$\psi(\tilde{x})/r_\psi$','interpreter','Latex','FontSize',labelsize)
ylim([0 1])
xlim([0,10])
ylim([0 1.05])
% plot(t(1:plotInd),U(1,1:plotInd),'k-')

% Error hist
subtightplot(3,2,5,subPlotGap,subPlotH,subPlotW);
if ws_flag
legend(h_phi,'$r_\phi$','interpreter','Latex','Fontsize',legendsize,'location','southwest')
end
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$\phi(\tilde{x},\hat{\textnormal{u}})$','interpreter','Latex','FontSize',labelsize)
% ylim([1e-12 1e-2])
if ws_flag
    ylim([1e-12 Inf])
    yticks([1e-12 1e-9 1e-6])
else
    yticks([1e-4 1e-2 1e0])
end
xlim([0,10])
ylim([1e-5 Inf])
yticks([1e-5 1e-4 1e-3 1e-2 1e-1])

% Iterations
subtightplot(3,2,6,subPlotGap,subPlotH,subPlotW);
xlabel('Time [s]','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell$','interpreter','Latex','FontSize',labelsize)
% yticks([0 10 20 30])
% legend(h_ell,'$\ell^*_{\textrm{ws}}$','interpreter','Latex','Fontsize',legendsize,'location','northeast')
xlim([0,10])
yticks([1 1e1 1e2 1e3 1e4])


if saveFigFlagOuter == 1
    filename = ['./Plots/traj_pend',num2str(ws_flag)];
    saveas(gcf,filename,'epsc')
end

%% Plot ROA business



% Plot ellStar trend
figure
load('./Data/ellStarTrend.mat')
figSize = [0 0 0.35 0.2]*.75;
set(gcf,'units','normalized','position',figSize)
% subPlotH = [0.22 0.04];
% subtightplot(1,2,2,subPlotGap,subPlotH,subPlotW)
semilogy(NVec,ellStar_cs,'k-')
hold on; box on; grid on;
semilogy(NVec_peak,ellPeak_cs,'k-.')
semilogy(NVec_peak,ellPeak_cs,'.','color',[0 0 0],'markersize',10)
 legend('$\ell^*$','$\textrm{max}(\ell_k)$','interpreter','Latex','Fontsize',legendsize,'location','southeast','NumColumns',1)
xlim([5 40])
% ylim([1 5e3])
xlabel('$N$','interpreter','Latex','Fontsize',labelsize)
ylabel('Iterations','interpreter','Latex','FontSize',labelsize)
% semilogy(NVec(1:iPlot),ellStar_cs(1:iPlot)-ellStar_ws(1:iPlot),'color',[0 0 0])
yticks([1 1e1 1e2 1e3 1e4 1e5])


if saveFigFlagOuter == 1
    filename = strcat('./Plots/NTrend_pend');
    saveas(gcf,filename,'epsc')
end



%}



