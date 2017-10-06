function ctd = getFlu(ctd,vvar);

%coeffs good for both pre and post April 2008
scale = 15.5; % ug/L/V
Voff = -0.053; % V
ctd.Flu = scale*(ctd.(vvar)-Voff);
