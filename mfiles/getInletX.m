function [x,y]=getInletX(lon,lat,doplot);

if nargin<3
  doplot=0;
end;

% define a Saanich Inlet co-ord system...
lon0=-123.5;
lat0=48.65;

% SJT - add mpernm to get code to work
mpernm = 1852;
%
x = (lon-lon0)*cos(lat0*pi/180)*mpernm*60+...
  sqrt(-1)*(lat-lat0)*mpernm*60;


xx =[
  -0.288709677419355  -1.504385964912280
  -0.321889400921659  -1.171052631578947
  -0.081336405529954  -0.960526315789473
   0.022350230414747   0.565789473684211
   0.781336405529954   1.021929824561404];
 xx = xx*1e4;
 xx = xx(:,1)+sqrt(-1)*xx(:,2);
 dx = 0*xx;
 dx(1)=0;
 
 for ii=1:length(xx)-1
   x0 = xx(ii);
%   rot = angle(diff(xx(ii+[0 1])));
   
 %  xnew(:,ii)=(x-x0)*exp(-sqrt(-1)*rot);
 %  dx(ii+1) = real((xx(ii+1)-x0)*exp(-sqrt(-1)*rot));
   dx(ii+1)=abs(xx(ii+1)-x0);
   dx(ii+1) = dx(ii+1)+dx(ii);
   
 end;
 xx = interp1(dx,xx,linspace(0,max(dx),1000));
 % smooth it
 N = 40;
 xn = conv2(xx,ones(1,N)/N,'same');
 xx((N+1):end-N) = xn((N+1):end-N);
 
 dx = cumsum([0 abs(diff(xx))]);
 
 x0=x;
 y = NaN*lon;
 x = NaN*lat;
 for i=1:length(lon);
   [yy,ind]=min(abs(xx-x0(i)))  
   
   % figure out if to the right or left.  right will be positive...
   if ind>=length(xx);
     ind = ind-1;
   end;
   an=angle(xx(ind+1)-xx(ind));     
   xn = (x0(i)-xx(ind))*exp(sqrt(-1)*-an);
   
   
   
   x(i) = real(xn)+dx(ind);
   y(i) = imag(xn);
   
 end;

 if doplot
 
   figure(1);clf
   load ../topo/SouthVanIsle
   xtopo=(VanIsleTopo.Lon-lon0)*cos(lat0*pi/180)*mpernm*60;
   ytopo=(VanIsleTopo.Lat-lat0)*mpernm*60;
   pcolor((xtopo)/1e3,(ytopo)/1e3,VanIsleTopo.z);
   shading flat;
   caxis([-1 1]*200);
   colormap([seajmk3(32);land3(32)]);
   hold on;
   plot(xx/1e3,'linewi',2);
   hold on;
   axis([-1 1 -1 1]*20)
   set(gca,'dataaspectRatio',[1 1 1])
   load ../200705/ctd/CtdGrid
   x=(cgrid.lon-lon0)*cos(lat0*pi/180)*mpernm*60;
   y=(cgrid.lat-lat0)*mpernm*60;
   
   plot(x/1e3,y/1e3,'d','markersize',6,'markerfacecol','k');
   
   for i=1:length(x);
     num = cgrid.id{i};
     text(x(i)/1e3,y(i)/1e3,sprintf(' %s\n x = %1.2f\n y = %1.2f',num,cgrid.xinlet(i)/1e3,cgrid.yinlet(i)/1e3),'fontsi',8);
   end
   xlabel('X [km]');
   ylabel('Y [km]');
   print -djpeg95 -r200 InletXMap.jpg
   unix('convert -geometry 700 InletXMap.jpg InletXMapSm.jpg');
   
 end
 
 