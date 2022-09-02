
close all; clear; clc

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

% Plot airfoil data (underformed & deformed)
figure; hold on;
plot(data(1:66,1),data(1:66,2),'-k','linewidth',1.5)
plot(data(67:end,1),data(67:end,2),'-k','linewidth',1.5)
plot(dataDef(1:66,1),dataDef(1:66,2),'-k','linewidth',1.5)
plot(dataDef(67:end,1),dataDef(67:end,2),'-k','linewidth',1.5)

box off
axis equal; %axis off
axis([-0.5 1.5 -1 0.5])


% Plot Elastic Axis & Hinge Axis Points
plot(0.532,yCenter,'ko','markerfacecolor','k','markersize',7)
plot(0.532,yDef,'ko','markerfacecolor','k','markersize',7)
elasticAxis = annotation('textarrow', [0.65 0.54], [0.80 0.655],'String','Eixo Elastico','interpreter','latex','FontSize',12);

plot(0.75, 0,'ko','markerfacecolor','k','markersize',7)
hingeAxis = annotation('textarrow', [0.70 0.625], [0.50 0.635],'String',{'Eixo de Rotacao da','Superficie de Controle'},'interpreter','latex','FontSize',12);


% Plot & Show variables
plot([1 1],[-0.05 0.05],'-.k','linewidth',0.1)
plot([0 0],[-0.05 0.05],'-.k','linewidth',0.1)
text(0,0.08,'$-b$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
text(1,0.08,'$b$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

plot([0.532 0.532],[-0.08 0.08],'-.k','linewidth',0.1)
text(0.532 ,0.12,'$ab$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

plot([0.75 0.75],[-0.06 0.06],'-.k','linewidth',0.1)
text(0.75 ,0.08, '$eb$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')


% Plot Reference Lines & Axis
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
text(0.40,  0.36, '$y$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

% Plot DOFs variables & geometric representation
pluge_Note = annotation('arrow', [0.53 0.53], [0.647 0.25]);
text(0.57, -0.50, '$h$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

theta_Note = annotation('arrow', [0.30 0.31], [0.235 0.31]);
text(-0.11, -0.73, '$\theta$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')


% Plot Undisturbed Velocity Arrows
VelocityArrow1 = annotation('arrow', [0.05 0.20], [0.648 0.648]); 
VelocityArrow2 = annotation('arrow', [0.05 0.20], [0.235 0.235]); 
text(-0.40,  0.08, '$V_{\infty}$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')
text(-0.40, -0.72, '$V_{\infty}$','interpreter','latex','FontSize',12,'HorizontalAlignment','center')

% Description
text(0.50, -1.1, 'Deformado','interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment','center')
text(0, 0.25, 'Indeformado','interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment','center')

axis off