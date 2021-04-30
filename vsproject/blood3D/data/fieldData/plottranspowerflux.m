plot(VarName4)
hold on
plot(VarName1,'--')
hold on
plot(VarName2,':')
hold off

xlabel('depth 0 - 12.1 um') % x-axis label
ylabel('Relative transmitted powerflux') % y-axis label
title('Transmitted powerflux - eps_rbc/eps_ba = 1.0993')
legend('epsi = 1.849 x 10^-7', 'epsi = 1.849 x 10^-9','epsi = 1.849 x 10^-11')
grid on

figure
plot(20*log10(VarName4))
hold on
plot(20*log10(VarName1),'--')
hold on
plot(20*log10(VarName2),':')
hold off
xlabel('depth 0 - 12.1 um') % x-axis label
ylabel('Relative transmitted powerflux (dB)') % y-axis label
title('Transmitted powerflux - eps_rbc/eps_ba = 1.0993')
legend('epsi = 1.849 x 10^-7', 'epsi = 1.849 x 10^-9','epsi = 1.849 x 10^-11')
grid on