function ctd2mat(ddir,coffset,exten);
% function ctd2mat(ddir);
% try to make a set of ctd files...


fnm = dir([ddir '/*.' exten])


for i=1:length(fnm);
  clear ctd

  if exten == 'hex'
    ctd = CtdHex2mat([ddir '/' fnm(i).name]);
    % find some useful tidbits if they occur
    ctd.lon=[];
    ctd.lat=[];
    ctd.id = '';
    ctd.time = [];
    for jj=1:length(ctd.headinfo)
      st = ctd.headinfo{jj}

      if length(st)>6 & strcmp(st(1:6),'** Lat');
        in = find(st==':');
        if ~isempty(in)
          st = st(in+1:end);
        else
          while ~ismember(st(1),char((0:9)+'0'))
            st=st(2:end);
          end;
        end;
        [ll]=sscanf(st,'%d %f%c');
        ctd.lat = ll(1)+ll(2)/60;
       end;

      if length(st)>6 & strcmp(st(1:6),'** Lon');
        printf('%s\n',st)
        printf('%s\n',st)
        [0:9]+'0'
        while ~ismember(st(1),char((0:9)+'0'))
          st=st(2:end);
        end;
        [ll]=sscanf(st,'%d %f%c');
        lon = ll(1)+ll(2)/60;
        ctd.lon = -lon;
      end
      if length(st)>6 & strcmp(st(1:6),'* Seac');
        [xx,dat] = sscanf(st,'* SeacatPlus V 1.6a  SERIAL NO. %d   %s');
        ctd.seriannum = xx(1);
        ctd.time =  datenum(st(end-20:end));
      end
      if length(st)>6 & strcmp(st(1:6),'** Sta');
        [id] = sscanf(st,'** Station ID:%s');
        ctd.id=id;
      end;
      ctd = getO2(ctd,'v1');
      ctd = getFlu(ctd,'v2');
      ctd = getPar(ctd,'v3');
    end;
  else
    clear ctd
    [lat,lon,gtime,data,names,sensors]=cnv2mat([ddir '/' fnm(i).name]);
    ctd.lon=lon;
    ctd.lat=lat;
    ctd.time = datenum(gtime);
    d = size(data)
    ctd.t = data(:,4);
    ctd.c = data(:,2);
    ctd.p = data(:,5);
    ctd.O2 = data(:,3);
  end;
  [dd,f]=fileparts(fnm(i).name);




  % find some useful tidbits if they occur

  % offset the c sensor...
  if nargin<2
    coffset=0;
  end;
  ctd.coffset=coffset;
  d = size(ctd.c)
  d = length(ctd.c)
  if ~isfield(ctd,'c0');
    ctd.c0 = ctd.c;
    ctd.c = interp1(1:length(ctd.c),ctd.c,(1:length(ctd.c))- ...
                    coffset)';
    d=length(ctd.c)
    ctd.sal = sw_salt(ctd.c*10/sw_c3515,ctd.t,ctd.p);
  end;

  [ddir '/' f '.mat']
  save('-7',[ddir '/' f '.mat'],'ctd')

end;
