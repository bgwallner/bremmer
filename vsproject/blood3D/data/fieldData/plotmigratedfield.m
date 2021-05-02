% pcolorplot of matrix
subplot(2,2,1)
s=pcolor(migratedfield16);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')

subplot(2,2,2)
s=pcolor(migratedfield43);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')

subplot(2,2,3)
s=pcolor(migratedfield85);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')

subplot(2,2,4)
s=pcolor(migratedfield128);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
