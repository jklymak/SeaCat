function seacat = seacatHex2mat(fname)
% function dat = seacatHex2mat(fname,offset)
% reads from a home2002 style seacat file.
%

  fin = fopen(fname,'r');
  
  if fin<0
    error(sprintf('Could not open %s',fname));
  end;

  % time
  l = fgetl(fin);
  hh = 0;
  while ~feof(fin) & (isempty(l) | l(1)=='*');
    hh = hh+1;
    seacat.headinfo{hh} = l;
    l
    try 
       i = strfind(l,'=');
       if ~isempty(i)
         v=l(2:i-1);
         val = l(i+1:end);
         eval(['seacat.head.CTDCoeff.' lower(v) '=' val ';']);
       end;
    end;
    l = fgetl(fin);
  end;
  
  t=fread(fin,Inf,'uchar');
  printf('%d',length(t))
  % t(end-3:end)

  printf('%d',length(t))
  

  s = reshape(t,36,[])';
  s(:,35:36)=[];  
  seacat.ctd = char(s);
  
  seacat = getTemp(seacat);
  seacat = getPres(seacat);
  seacat = getCond(seacat);
  seacat = getVolt(seacat);
  
  %keyboard;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getTemp;
function fctd = getTemp(fctd);

  a0 = fctd.head.CTDCoeff.ta0;
  a1 = fctd.head.CTDCoeff.ta1;
  a2 = fctd.head.CTDCoeff.ta2;
  a3 = fctd.head.CTDCoeff.ta3;
  
  rawT = hex2dec(fctd.ctd(:,1:6));
  mv = (rawT-524288)/1.6e7;
  r = (mv*2.9e9 + 1.024e8)./(2.048e4-mv*2e5);
  fctd.t = a0+a1*log(r)+a2*log(r).^2+a3*log(r).^3;
  fctd.t = 1./fctd.t - 273.15;
  return;

function fctd = getCond(fctd);

  g = fctd.head.CTDCoeff.g;
  h = fctd.head.CTDCoeff.h;
  i = fctd.head.CTDCoeff.i;
  j = fctd.head.CTDCoeff.j;
  tcor = fctd.head.CTDCoeff.ctcor;
  pcor = fctd.head.CTDCoeff.cpcor;
  
  f = hex2dec(fctd.ctd(:,(1:6)+6))/256/1000;
  
  fctd.c = (g+h*f.^2+i*f.^3+j*f.^4)./(1+tcor*fctd.t+pcor*fctd.p);
  
  return;
function fctd = getVolt(fctd);

  fctd.v1 = hex2dec(fctd.ctd(:,(23:26)))/13107;
  fctd.v2 = hex2dec(fctd.ctd(:,(23:26)+4))/13107;
  fctd.v3 = hex2dec(fctd.ctd(:,(23:26)+8))/13107;
  
  
  return;

function fctd = getPres(fctd);

  pa0 = fctd.head.CTDCoeff.pa0;
  pa1 = fctd.head.CTDCoeff.pa1;
  pa2 = fctd.head.CTDCoeff.pa2;
  ptempa0 = fctd.head.CTDCoeff.ptempa0;
  ptempa1 = fctd.head.CTDCoeff.ptempa1;
  ptempa2 = fctd.head.CTDCoeff.ptempa2;
  ptca0 = fctd.head.CTDCoeff.ptca0;
  ptca1 = fctd.head.CTDCoeff.ptca1;
  ptca2 = fctd.head.CTDCoeff.ptca2;
  ptcb0 = fctd.head.CTDCoeff.ptcb0;
  ptcb1 = fctd.head.CTDCoeff.ptcb1;
  ptcb2 = fctd.head.CTDCoeff.ptcb2;
  
  rawP = hex2dec(fctd.ctd(:,(1:6)+12));

  y = hex2dec(fctd.ctd(:,(1:4)+18))/13107;
  %  y = 524965.58
  t = ptempa0+ptempa1*y+ptempa2*y.^2;
  x = rawP-ptca0-ptca1*t-ptca2*t.^2;
  n = x*ptcb0./(ptcb0+ptcb1*t+ptcb2*t.^2);
  
  fctd.p = (pa0+pa1*n+pa2*n.^2-14.7)*0.689476;