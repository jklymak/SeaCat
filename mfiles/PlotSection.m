load data/CtdGrid
M = size(cgrid.c,2)
cgrid.sal=sw_salt(cgrid.c*10/sw_c3515,cgrid.t,repmat(cgrid.depths,1,M))
cgrid.pden = sw_pden(cgrid.sal,cgrid.t,repmat(cgrid.depths,1,M),0)


%save data/CtdGridW cgrid




%%
clf

t0=datenum(2014,1,1)
subplot2(4,1,1)
facetsurface(cgrid.alongx,cgrid.depths,cgrid.t);shading flat
caxis([8,15])
hold on
contour(cgrid.alongx,cgrid.depths,cgrid.pden-1000.,17:25,'k')
axis ij;
colorbar()
title('T [^oC]')
%xlim([188.65,189])
subplot2(4,1,2)
facetsurface(cgrid.alongx,cgrid.depths,cgrid.sal);shading flat
caxis([20,32])
hold on
contour(cgrid.alongx,cgrid.depths,cgrid.pden-1000.,17:26,'k')
axis ij;
colorbar()
title('S [psu]')
%xlim([188.65,189])
subplot2(4,1,3)
facetsurface(cgrid.alongx,cgrid.depths,cgrid.O2);shading flat
caxis([0,7])
hold on
contour(cgrid.alongx,cgrid.depths,cgrid.pden-1000.,17:0.5:28,'k')
axis ij;
colorbar()
title('O2 [mL/L]')
%xlim([188.65,189])
subplot2(4,1,4)
facetsurface(cgrid.alongx,cgrid.depths,cgrid.Flu);shading flat
%caxis([0,7])
hold on
contour(cgrid.alongx,cgrid.depths,cgrid.pden-1000.,17:25,'k')
axis ij;
colorbar()
title('Flu')

%% Plot T/S
jmkfigure(2,2,0.4)
clf
subplot(1,2,1)
plot(cgrid.t,cgrid.sal,'.')
xlim([8,15])
ylim([27,32.5])
axis ij

subplot(1,2,2)
plot(cgrid.O2,cgrid.sal,'.')
%xlim([8,15])
ylim([27,32.5])
axis ij