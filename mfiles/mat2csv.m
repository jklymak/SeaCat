function mat2csv(ddir);

d = dir([ddir '/*.mat']);

for i=1:length(d);
  nm = d(i).name(1:end-4);
  
  clear ctd;
  load([ddir d(i).name]);
  if exist('ctd','var');
    nm = d(i).name(1:end-4)
    ctd.sal = sw_salt(ctd.c*10/sw_c3515,ctd.t,ctd.p);
    ctd.pden=sw_pden(ctd.sal,ctd.t,ctd.p,0)-1000;
    M = [ctd.p ctd.t ctd.c ctd.sal ctd.pden ctd.O2 ctd.O2new ctd.Flu log10(ctd.Par) ctd.v1 ctd.v2 ctd.v3];
    fin = fopen([ddir '/' nm '.csv'],'w');


    fprintf(fin,'fname, %s',d(i).name);
    fprintf(fin,'\n');

    fprintf(fin,'Lon, %f',ctd.lon);
    fprintf(fin,'\n');

    fprintf(fin,'Lat, %f',ctd.lat);
    fprintf(fin,'\n');
    fprintf(fin,'Time, %s\n',datestr(ctd.time));

    if isfield(ctd,'id');
      fprintf(fin,'Id, %s',ctd.id);
    end;
    fprintf(fin,'\n');

    fprintf(fin,'P [dbar] ,T [^oC],C [S/m],S [psu],PotDens [kg/m^3],02 [mg/L], O2 [umol/kg], Flouresence [uG/L],log10(Par [uEinsteins/m^2/s]),V1 [V],V2 [V],V3 [V]\n');
    fclose(fin);
    dlmwrite([ddir '/' nm '.csv'],M,'-append','precision','%6.6f');
  end;
end;