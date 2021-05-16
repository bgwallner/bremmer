% pcolorplot of matrix

figure
s=pcolor(migratedfield32);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
title('Absolute electric field in transverse direction. Euler angle \theta=0. z=1.9\mum')
xlabel('Transverse - 0 < x < 62 \mum')
ylabel('Transverse - 0 < y < 62 \mum') 

figure
s=pcolor(migratedfield64);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
title('Absolute electric field in transverse direction. Euler angle \theta=0. z=3.9\mum')
xlabel('Transverse - 0 < x < 62 \mum')
ylabel('Transverse - 0 < y < 62 \mum') 

figure
s=pcolor(migratedfield96);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
title('Absolute electric field in transverse direction. Euler angle \theta=0. z=5.8\mum')
xlabel('Transverse - 0 < x < 62 \mum')
ylabel('Transverse - 0 < y < 62 \mum')

figure
s=pcolor(migratedfield128);
s.FaceColor='interp';
set(s, 'EdgeColor', 'none')
colorbar
caxis([0 1.5])
title('Absolute electric field in transverse direction. Euler angle \theta=0. z=7.8\mum')
xlabel('Transverse - 0 < x < 62 \mum')
ylabel('Transverse - 0 < y < 62 \mum')

