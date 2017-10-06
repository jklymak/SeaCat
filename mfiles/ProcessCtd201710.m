addpath /Users/jklymak/mixingsoftware/seawater
addpath /Users/jklymak/mixingsoftware/general
% the directory with the data in it
bdir = './201710/';

coffset = 1.7; % conductivity offset in scans.

ctd2mat(bdir,coffset, 'hex');
ls 201709

% any fixing gets done here

% make a csv file of the mat file
% mat2csv(bdir);

% bin the ctd data into 1-m bins.
binCtd([bdir])
load([bdir 'CtdGrid.mat'])

%%%
% Re-arrange any casts that are out of order
if 0
	ncasts=size(cgrid.t,2)
	ind = 1:ncasts

	ind(13:16)=[16 13 14 15]
	% ind = fliplr(ind)

	fnms=fieldnames(cgrid)
	for ii = 1:length(fnms)
	  if size(cgrid.(fnms{ii}),2)==ncasts
	    fnms{ii}
	    cgrid.(fnms{ii})=cgrid.(fnms{ii})(:,ind);
	  end
	end
end
%% get along-track from S4

cgrid=rmfield(cgrid,'xinlet')
cgrid=rmfield(cgrid,'yinlet')

cgrid.alongx(1)=0.;
for i=1:length(cgrid.lon)-1

  cgrid.alongx(i+1)=cgrid.alongx(i)+sw_dist(cgrid.lat(i+[0,1]),cgrid.lon(i+[0,1]),'km');
end

% OK, x is backwards
cgrid.alongx=cgrid.alongx(end)-cgrid.alongx;

save('-7',[bdir,'CtdGrid.mat'], 'cgrid')
