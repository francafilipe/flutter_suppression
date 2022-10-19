clear; close all; clc;

%% 
k = [0:0.01:10]; 
Ck = theodorsen(k);

figure
subplot(2,1,1); hold on; box on;
plot(k,real(Ck),'-k','linewidth',2)
xlabel('$k = \frac{\omega b}{V_{\infty}}$ [-]','interpreter','latex','FontSize',12); ylabel('amplitude de $C(k)$ [-]','interpreter','latex','FontSize',12);
subplot(2,1,2); hold on; box on;
plot(k,phase(Ck)*180/pi,'-k','linewidth',2)
xlabel('$k = \frac{\omega b}{V_{\infty}}$ [-]','interpreter','latex','FontSize',12); ylabel('fase de $C(k)$ [$^{\circ}$]','interpreter','latex','FontSize',12);
