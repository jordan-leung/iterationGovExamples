clc
clear all
close all

set(0,'defaultLineLineWidth',1)
set(0,'defaultAxesFontSize',12)


%% Plot trajectories (CS)
saveFigFlagOuter = 0;
labelsize = 14;
legendsize = 12;

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

figure
set(gcf,'units','normalized','position',figSize)
icflag = 3;
dataList =[{strcat('output_ig_N10_ic',num2str(icflag),'_ws0')};...
    {strcat('output_sg_N10_ic',num2str(icflag),'_ws0')};...
    {strcat('output_tracking_N46_ic',num2str(icflag))}];
handVec = zeros(3,1);
handVec2 = zeros(3,1);

% s
for iOuter = 1:3
    % Load IG
    loadstr = ['./Data/',dataList{iOuter}];
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
    handVec(iOuter) = plot(t,X(1,:),'Color',colorMatrix(iOuter,:));
    hold on; box on; grid on
    if iOuter < 3
        plot(t,refHist,'Linestyle','-.','LineWidth',1,'Color',colorMatrix(iOuter,:))
    end
%     plot(t(1),X(1,1),markerStr{iOuter},'MarkerSize',sillyMarkerSize,'Color',colorMatrix(iOuter,:),'linewidth',1);
%     plot(t(1),X(1,1),'.','Color',colorMatrix(iOuter,:),'linewidth',1);

    % delta
    subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
    if iOuter == 1
        h_u = plot(t,t*0+umax(1),'k-.','LineWidth',1);
        hold on; box on; grid on;
        plot(t,t*0+umin(1),'k-.','LineWidth',1)
    end
    plot(t(1:end-1),U(1,:),'LineWidth',1,'Color',colorMatrix(iOuter,:))


    % Iterations (CG and IG)
    if iOuter < 3
        subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW);
        if iOuter == 1
            h_ell = plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.ellStar_cs,'k-.');
            hold on; box on; grid on;
        end
        handVec2(iOuter) = stairs(t(1:end-1),ellHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
    else
        subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
        h_ell2 = plot(t(1:end-1),t(1:end-1)*0+output.controlArgs.const.ellStar_cs,'k-.');
        hold on; box on; grid on;
        handVec(iOuter) = stairs(t(1:end-1),ellHist,'LineWidth',1,'color',colorMatrix(iOuter,:));
    end 
end
% Add legends and axes
% x1
subtightplot(3,2,1,subPlotGap,subPlotH,subPlotW)
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$s$','interpreter','Latex','FontSize',labelsize)
% legend('$x_{1,k}$','$v_k$','interpreter','Latex','Fontsize',legendsize,'location','northeast')
ylim([-1.1 1.1])
legend(handVec,'IG-MPC','SG-MPC','rMPC','interpreter','Latex','Fontsize',legendsize,'location','southeast')


% u 
subtightplot(3,2,2,subPlotGap,subPlotH,subPlotW)
% legend(h_u,'$\mathcal{U}$','interpreter','Latex','Fontsize',legendsize,'location','southeast')
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$u$','interpreter','Latex','FontSize',labelsize)
ylim([-0.105 0.105])

% Iterations
subtightplot(3,2,3,subPlotGap,subPlotH,subPlotW);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell$','interpreter','Latex','FontSize',labelsize)
% ylim([0 35])
% yticks([0 10 20 30])
legend([handVec2(1:2); h_ell],'IG-MPC','SG-MPC','$\ell^*_{N=10}$','interpreter','Latex','Fontsize',legendsize,'location','northeast')

% Iterations
subtightplot(3,2,4,subPlotGap,subPlotH,subPlotW);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\ell$','interpreter','Latex','FontSize',labelsize)
% ylim([0 35])
% yticks([0 10 20 30])
legend([handVec(3); h_ell2],'rMPC','$\ell^*_{N=46}$','interpreter','Latex','Fontsize',legendsize,'location','northeast')

if saveFigFlagOuter == 1
    figure(1)
    filename = strcat('./Plots/IGvSG-Compare_doubleInt');
    saveas(gcf,filename,'epsc')
end
