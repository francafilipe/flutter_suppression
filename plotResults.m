function plotResults(T,k,x,input)
%% 

%%
% Define time array to plot
t = linspace(0,T,k+1);

% State Vector Dynamics
figure(1)

subplot(2,2,1);
plot(t,x(1,:),'-k','linewidth',1.5)
ylabel('$h [in]$','Interpreter','latex'); 
xlabel('$time [sec]$','Interpreter','latex');

subplot(2,2,2);
plot(t,x(3,:),'-k','linewidth',1.5)
ylabel('$\dot{h} [in/s]$','Interpreter','latex'); 
xlabel('$time [sec]$','Interpreter','latex');

subplot(2,2,3);
plot(t,x(2,:)*180/pi,'-k','linewidth',1.5)
ylabel('$\alpha [deg]$','Interpreter','latex'); 
xlabel('$time [sec]$','Interpreter','latex');

subplot(2,2,4);
plot(t,x(4,:)*180/pi,'-k','linewidth',1.5)
ylabel('$\dot{\alpha} [deg/s]$','Interpreter','latex'); 
xlabel('$time [sec]$','Interpreter','latex');


% TE Actuator & Control Input
figure(2)

subplot(2,1,1)
plot(t, input*180/pi, '-k', 'linewidth', 1.5)
ylabel('$u_{c,TEA} [deg]$','Interpreter','latex'); xlabel('$time [sec]$','Interpreter','latex');

subplot(2,1,2)
hold on;
plot(t,x(13,:)*180/pi,'-k','linewidth',1.5)
plot(t,x(14,:)*180/pi*1e-1,'--r','linewidth',1.5)
plot(t,x(15,:)*180/pi*1e-3,'-.b','linewidth',1.)
ylabel('$u_{TEA} [deg]$','Interpreter','latex'); xlabel('$time [sec]$','Interpreter','latex');
legend('$u_{TEA}$','$\dot{u}_{TEA}*10^{-1}$','$\ddot{u}_{TEA}*10^{-3}$', 'Interpreter', 'latex');


end