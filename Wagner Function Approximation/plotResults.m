function plotResults(T,k,x,input)
%% 

%%
% Define time array to plot
t = linspace(0,T,k+1);

% State Vector Dynamics
figure(1)

subplot(3,2,1);
plot(t,x(1,:),'-k','linewidth',1.5)
ylabel('h [m]'); xlabel('time [sec]');

subplot(3,2,2);
plot(t,x(4,:),'-k','linewidth',1.5)
ylabel('h_{dot} [m/s]'); xlabel('time [sec]');

subplot(3,2,3);
plot(t,x(2,:),'-k','linewidth',1.5)
ylabel('\alpha [rad]'); xlabel('time [sec]');

subplot(3,2,4);
plot(t,x(5,:),'-k','linewidth',1.5)
ylabel('\alpha_{dot} [rad/s]'); xlabel('time [sec]');

subplot(3,2,5);
plot(t,x(3,:),'-k','linewidth',1.5)
ylabel('\beta [rad]'); xlabel('time [sec]');

subplot(3,2,6);
plot(t,x(6,:),'-k','linewidth',1.5)
ylabel('\beta_{dot} [m]'); xlabel('time [sec]');


% Control Action
figure(2)
plot(t, input, '-k', 'linewidth', 1.5)
ylabel('u'); xlabel('time [sec]');



end