function ctd = getPar(ctd,vvar);

%factors good for both pre and post April 2008
ConFac = 10111223458.03842;
offset = -0.15423974;

ctd.Par =  (ctd.(vvar)-offset)*ConFac;