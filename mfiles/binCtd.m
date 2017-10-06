function binCtd(ddir);


d = dir([ddir '*.mat']);

%
clear cgrid;
num = 0;
for i=1:length(d);
  clear ctd;
  if ~strcmp(d(i).name,'CtdGrid.mat')

    load([ddir '/' d(i).name]);
    exist('ctd','var')
    if exist('ctd','var');

      in = find(diff(ctd.p)>0.02 & ctd.p(2:end)>4);
      zbin =[0:325];
      % ctd.sal = sw_salt(ctd.c*10/sw_c3515,ctd.t,ctd.p);

      if ~isfield(ctd,'c0');
        ctd.c0 = ctd.c;
        ctd.c = interp1(1:length(ctd.c),ctd.c,(1:length(ctd.c))-1.7)';
        save([ddir '/' d(i).name],'ctd');
      end;

      todo = {'t','c','O2'};% ,'Flu','Par'};

      if 1 % ~(ctd.lon<-123.52 & ctd.lat>48.72)
        num = num+1
        if (isfield(ctd,'lat') & length(ctd.lat)>0)
          cgrid.lat(num)=ctd.lat;
          cgrid.lon(num)=ctd.lon;
        else
          cgrid.lat(num)=NaN;
          cgrid.lon(num)=NaN;
        end;
        cgrid.name{num}=d(num).name;
        if isfield(ctd,'id');
          cgrid.id{num} = ctd.id;
        else
          cgrid.id{num}='';
        end
        cgrid.time(num)=ctd.time
  
        for j=1:length(todo);
          todo{j}
          cgrid.(todo{j})(:,num) = bindata1d(zbin,ctd.p(in),ctd.(todo{j})(in))';
          cgrid.depths = zbin(1:end-1)'+0.5*(zbin(2)-zbin(1));
        end;
      end;
    end
  end

end;


% get x....
try
  [cgrid.xinlet,cgrid.yinlet] = getInletX(cgrid.lon,cgrid.lat);
catch
  
end

[ll,ind]=sort(abs(cgrid.time));

fnames = fieldnames(cgrid);
for i=1:length(fnames);
  if ~strcmp(fnames{i},'name') & ~strcmp(fnames{i},'id')
    if size(cgrid.(fnames{i}),2)==length(ind);
      cgrid.(fnames{i})=cgrid.(fnames{i})(:,ind);
    end;
  end;
end;

nm = cgrid.name;
for i=1:length(ind);
  cgrid.name{i} = nm{ind(i)};
end;
nm = cgrid.id;
for i=1:length(ind);
  cgrid.id{i} = nm{ind(i)};
end;

save([ddir 'CtdGrid.mat'],'cgrid');

% OK make a csv file for each cast
% makeCsv(ddir,cgrid);

%%
function makeCsv(ddir,cgrid)

for i=1:length(cgrid.lon);
  nm = [cgrid.name{i}(1:end-4) '_1m.csv'];
  fin = fopen([ddir nm],'w');
  fprintf(fin,'fname, %s',nm(1:end-7));
  fprintf(fin,'\n');

  fprintf(fin,'Lon, %f',cgrid.lon(i));
  fprintf(fin,'\n');

  fprintf(fin,'Lat, %f',cgrid.lat(i));
  fprintf(fin,'\n');
  fprintf(fin,'Time, %s\n',datestr(cgrid.time(i)));
  fprintf(fin,'P [dbar] ,T [^oC],C [S/m],S [psu],PotDens [kg/m^3],O2 [mg/L], Fluoresence [uG/L],log10(Par [uEinsteins/m^2/s]),V1 [V],V2 [V],V3 [V]\n');
  fclose(fin);
  cgrid.sal=sw_salt(cgrid.c*10/sw_c3515,cgrid.t,repmat(cgrid.depths,1,size(cgrid.c,2)));
  cgrid.pden=sw_pden(cgrid.sal,cgrid.t,repmat(cgrid.depths,1,size(cgrid.t,2)),0);
  
  M = [cgrid.depths];
%  todo = {'t','c','sal','pden','O2','O2new','Flu','Par'};
  todo = {'t','c','sal','pden','O2','Flu','Par'};
for ii=1:length(todo);
    M = [M cgrid.(todo{ii})(:,i)];
    
  end
  dlmwrite([ddir '/' nm ],M,'-append','precision','%6.6f');

end;
