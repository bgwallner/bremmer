figure
plot(migratedfield580(:,512),'-')
set(gca, 'YScale', 'log')
xlim([512 750])
grid on

hold on
plot(migratedfield580(512,:), '-.')
set(gca, 'YScale', 'log')
xlim([512 750])

% theta = 0
plot(migratedfield3(:,512),'--')
set(gca, 'YScale', 'log')
xlim([512 750])
