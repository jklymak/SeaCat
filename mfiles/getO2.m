function ctd = getO2(ctd,vvar);

ctd.sal = sw_salt(ctd.c*10/sw_c3515,ctd.t,ctd.p);

ctd.O2sat = real(sw_satO2(ctd.sal,ctd.t));
%from April 2008 calibration
Voff = -0.4877;
Soc = 0.4401;
Boc = 0.0000;
Tcor = 0.0005;
Pcor = 1.35e-4;

%below numbers are true for CTD casts prior to April 2008
%Voff = -0.4893;
%Soc = 0.4227;
%Boc = 0.0000;
%Tcor = 0.001;
%Pcor = 1.35e-4;

ctd.O2 = Soc*(ctd.(vvar)+Voff).*exp(Tcor.*ctd.t).*ctd.O2sat.*exp(ctd.p*Pcor);
%  ctd = rmfield(ctd,'sal');
ctd = rmfield(ctd,'O2sat');

%sjt added code below
ctd.den=sw_dens(ctd.sal, ctd.t, ctd.p);
ctd.O2= 1000000*ctd.O2./(ctd.den*22.3916);
