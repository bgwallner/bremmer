figure
plot(transpowerfluxzero)
hold on
plot(transpowerfluxpifourth)
hold on
plot(transpowerfluxpihalf,'--')
hold on
plot(transpowerfluxrandom,':')
hold on

volumeFraction = 0.2134;
lambda = 10.437938144329898;
nRBC = 4.3*10^-6;
z=4*pi*(volumeFraction*nRBC/lambda)*[1:1:1020];
beer = exp(-z);

plot(beer,':')
hold off

xlabel('Propagation depth 0 - 62\mum')   % x-axis label
ylabel('Relative transmitted intensity') % y-axis label
title('Transmitted intensity - \epsilon_{rbc}/\epsilon_{ba} = 1.0993 (\lambda=632.8nm)')
legend('\theta=0','\theta=\pi/4', '\theta=\pi/2','\theta=random', 'I(z)=e^{(-4\pi*V_{rbc}*Im(n_{rbc}) / (V_{tot}\lambda)}')
grid on
