close all

% set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesFontSize',12)
labelsize = 15;

saveDirectory = 'C:\Users\Jordan Leung\Documents\MATLAB\PhD\19-08-26_HCWToyProblem';


colorMatrix = [0 0 0; 228 26 28; 55 126 184; 77 175 74]./255;
greyColor = [111 111 111]/255;

%% Producing plots

% Unpack variables
X = output.X;
U = output.U;
refHist = output.refHist;
VHist = output.VHist;
ellHist = output.ellHist;

% Find corresponding regular time index
figSize = [0 0 0.3 0.3];
subPlotGap = [0.13 0.10];
subPlotH = [0.1 0.05];
subPlotW = [0.1 0.05];


figure
% s
set(gcf,'units','normalized','position',figSize)
subtightplot(2,3,1,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on
plot(t,X(1,:),'Color',colorMatrix(1,:))
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('${s}$','interpreter','Latex','FontSize',labelsize)
plot(t,r,'Linestyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(t,refHist,'Linestyle','-.','LineWidth',2,'Color',colorMatrix(3,:))
legend('${s}$','$r$','$v$','interpreter','Latex','Fontsize',12,'location','southeast')

% \beta
subtightplot(2,3,2,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,X(2,:),'Color',colorMatrix(1,:))
% plot(t,t*0+xmax,'-.','color',greyColor);
% plot(t,t*0+xmin,'-.','color',greyColor);
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\dot{s}$','interpreter','Latex','FontSize',labelsize)

% delta
subtightplot(2,3,3,subPlotGap,subPlotH,subPlotW)
hold on; box on; grid on;
plot(t,t*0+umax(1),'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t,t*0+umin(1),'Linestyle','-.','LineWidth',2,'Color',greyColor)
plot(t(1:end-1),U(1,:),'k','LineWidth',2)
xlim([0 max(t)])
xlabel('Time','interpreter','Latex','Fontsize',labelsize)
ylabel('$\delta$','interpreter','Latex','FontSize',labelsize)

% V
subtightplot(2,3,4,subPlotGap,subPlotH,subPlotW);
plot(t(1:end-1),VHist,'k','LineWidth',2)
hold on; box on; grid on;
plot(t,t*0+VBar,'r--','linewidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','Latex','Fontsize',labelsize)
ylabel('$V_N^*$','interpreter','Latex','FontSize',labelsize)
legend('$V(x,v)$','$N\gamma + \alpha$','interpreter','Latex','Fontsize',10,'location','best')


% ell
subtightplot(2,3,5,subPlotGap,subPlotH,subPlotW);
semilogy(t(1:end-1),ellHist,'k','LineWidth',2)
hold on; box on; grid on;
plot(t,t*0+ellStar_cs,'r--','linewidth',2)
xlim([0 max(t)])
xlabel('$t$','interpreter','Latex','Fontsize',labelsize)
ylabel('$V_N^*$','interpreter','Latex','FontSize',labelsize)
legend('$V(x,v)$','$N\gamma + \alpha$','interpreter','Latex','Fontsize',10,'location','best')





% if saveFigFlag == 1
%     filename = strcat('./Plots/output3');
%     saveas(gcf,filename,'epsc')
% end
% 
% 


