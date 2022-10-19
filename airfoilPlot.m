close all; clear; clc

%% PLOT AIRFOIL - UNDEFORMED & DEFORMED

% Airfoil Undeformed data points
data = importdata('airfoil.txt');

xCenter = 0.532;
yCenter = 0;
center = repmat([xCenter;yCenter],1,length(data(:,1)));

% Define Deformed data points
theta = -15;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

dataRot = (R*(data'-center) + center)';

yDef = -0.8;
centerDef = repmat([0;yDef],1,length(dataRot(:,1)));
dataDef = dataRot + centerDef';

% Define deflected flap points
flapData = dataDef(149:244,:);

xHinge = 0.79;
yHinge = -0.89;
hinge = repmat([xHinge;yHinge],1,length(flapData(:,1)));

delta = -25;
R_2 = [cosd(delta) -sind(delta); sind(delta) cosd(delta)];
flapRot = (R_2*(flapData'-hinge) + hinge)';

% Plot Airfoil
figure; hold on;
% Undeformed
plot(data(1:46,1),data(1:46,2),'-k','linewidth',1.5)
plot(data(47:92,1),data(47:92,2),'-k','linewidth',1.5)
plot(data(93:148,1),data(93:148,2),'-k','linewidth',1.5)
plot(data(149:168,1),data(149:168,2),'-k','linewidth',1.5)
plot(data(169:188,1),data(169:188,2),'-k','linewidth',1.5)
plot(data(189:244,1),data(189:244,2),'-k','linewidth',1.5)
% Deformed
plot(dataDef(1:46,1),dataDef(1:46,2),'-k','linewidth',1.5)
plot(dataDef(47:92,1),dataDef(47:92,2),'-k','linewidth',1.5)
plot(dataDef(93:148,1),dataDef(93:148,2),'-k','linewidth',1.5)
% plot(dataDef(149:168,1),dataDef(149:168,2),'--k','linewidth',.5)
% plot(dataDef(169:188,1),dataDef(169:188,2),'--k','linewidth',.5)
% plot(dataDef(189:244,1),dataDef(189:244,2),'--k','linewidth',.5)
% Flap Deflected
plot(flapRot(1:20,1),flapRot(1:20,2),'-k','linewidth',1.5)
plot(flapRot(21:40,1),flapRot(21:40,2),'-k','linewidth',1.5)
plot(flapRot(41:96,1),flapRot(41:96,2),'-k','linewidth',1.5)

% Plot Configuration
box off
axis equal; axis off
axis([-0.5 1.5 -1 0.5])


%% Reference Lines & Axis
horizontalChord = [[-0.15 1.15]; [0 0]];
plot(horizontalChord(1,:),horizontalChord(2,:),'--k')

horizontalChordDef = [[-0.15 1.15]; [yDef yDef]*1.00];
center = repmat([0.532;yDef],1,length(horizontalChordDef(1,:)));
deformedChord = R*(horizontalChordDef-center)+center;

plot(horizontalChordDef(1,:),horizontalChordDef(2,:),'--k')
plot(deformedChord(1,:),deformedChord(2,:),'--k')

xAxis_Note = annotation('arrow', [0.50 0.80], [0.647 0.647]); xAxis_Note.LineStyle = '-';
yAxis_Note = annotation('arrow', [0.50 0.50], [0.647 0.85]); yAxis_Note.LineStyle = '-';

text(1.20, -0.07, '$x$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
text(0.40,  0.36, '$z$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')


%% ELASTIC AXIS & ANNOTATIONS
plot(0.532,yCenter,'ko','markerfacecolor','k','markersize',7)
plot(0.532,yDef,'ko','markerfacecolor','k','markersize',7)
elasticAxis = annotation('textarrow', [0.65 0.54], [0.80 0.655],'String','Eixo El\''astico','interpreter','latex','FontSize',12);

plot([0.532 0.532],[-0.08 0.08],'-.k','linewidth',0.1)
text(0.532 ,0.12,'$ab$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')


%% HINGE AXIS & ANNOTATIONS
% plot(0.81, 0,'ko','markerfacecolor','k','markersize',7)
hingeAxis = annotation('textarrow', [0.70 0.64], [0.50 0.63],'String',{'Eixo de Rota\c{c}$\tilde{\textrm{a}}$o da','Superf\''icie de Controle'},'interpreter','latex','FontSize',12);

%% GDLs & DISPLACEMENTS
% Plunge (h) annotation
pluge_Note = annotation('arrow', [0.53 0.53], [0.647 0.25]);
text(0.48, -0.50, '$h$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

% Theta annotation
thetaArrowPoints = [-0.12 -0.80; -0.12 -0.718; -0.08 -0.635];
thetaArrowYPoints = thetaArrowPoints(1,2):0.01:thetaArrowPoints(3,2);
thetaArrowXPoints = spline(thetaArrowPoints(:,2),thetaArrowPoints(:,1),thetaArrowYPoints);
plot(thetaArrowXPoints,thetaArrowYPoints,'-k','linewidth',0.1)
text(-0.16, -0.73, '$\theta$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
theta_Note = annotation('arrow', [0.287 0.29], [0.31 0.32]);

% Beta annotation
betaArrowPoints = [0.9394 -0.9836; 0.955 -0.975; 0.98 -0.92];
betaArrowYPoints = betaArrowPoints(1,2):0.001:betaArrowPoints(3,2);
betaArrowXPoints = spline(betaArrowPoints(:,2),betaArrowPoints(:,1),betaArrowYPoints);
plot(betaArrowXPoints,betaArrowYPoints,'-k','linewidth',0.1)
text(1.01, -1.00, '$\beta$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
beta_Note = annotation('arrow', [0.698 0.69], [0.15 0.14]);

%% VARIABLES
% Plot & Show variables
plot([1 1],[-0.05 0.05],'-.k','linewidth',0.1)
plot([0 0],[-0.05 0.05],'-.k','linewidth',0.1)
text(0,0.08,'$-b$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
text(1,0.08,'$b$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

plot([0.532 0.532],[-0.08 0.08],'-.k','linewidth',0.1)
text(0.532 ,0.12,'$ab$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

plot([0.80 0.80],[-0.06 0.06],'-.k','linewidth',0.1)
text(0.75 ,0.08, '$eb$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

%% UNDISTURBED VELOCITY
VelocityArrow1 = annotation('arrow', [0.05 0.20], [0.648 0.648]); 
VelocityArrow2 = annotation('arrow', [0.05 0.20], [0.235 0.235]); 
text(-0.40,  0.08, '$V_{\infty}$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
text(-0.40, -0.72, '$V_{\infty}$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

%% DESCRIPTION
text(0.50, -1.1, 'Deformado','interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment','center')
text(0, 0.25, 'N$\tilde{\textrm{a}}$o\ Deformado','interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment','center')
