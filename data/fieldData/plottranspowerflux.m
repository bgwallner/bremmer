plot(VarName4)
hold on
plot(VarName1)
hold on
plot(VarName2,'--')
hold on
plot(VarName5,':')
hold off

xlabel('depth 0 - 12.1 um') % x-axis label
ylabel('Relative transmitted powerflux') % y-axis label
title('Transmitted powerflux - \epsilon_{rbc}/\epsilon_{ba} = 1.0993 (\lambda=632.8nm, random RBC)')
legend('epsi = 1.849 x 10^-11','epsi = 1.849 x 10^-7', 'epsi = 1.849 x 10^-5','epsi = 1.849 x 10^-3')
grid on

figure
plot(20*log10(VarName4))
hold on
plot(20*log10(VarName1))
hold on
plot(20*log10(VarName2),'--')
hold on
plot(20*log10(VarName5),':')
hold off
xlabel('depth 0 - 12.1 um') % x-axis label
ylabel('Relative transmitted powerflux (dB)') % y-axis label
title('Transmitted powerflux - eps_rbc/eps_ba = 1.0993')
legend('epsi = 1.849 x 10^-11','epsi = 1.849 x 10^-7', 'epsi = 1.849 x 10^-5','epsi = 1.849 x 10^-3')
grid on