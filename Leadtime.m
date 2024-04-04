%%%%%%%%%% Lead time to landfall %%%%%%%%%%%%
 
clear all;
close all; m_proj('equidist','lat',[20 50],'lon',[-99 -60],'sphere','aspect',10);

#### bathymetry data to calculate landfall #######

file = 'ETOPO_2022_v1_60s_N90W180_bed.nc'
ncdisp(file)

SID={'canesm','cnrm6','ecearth6','ipsl6','miroc6','ukmo6'}
LEADTIME_A=[];

for v = 1:6; % example models 
    op = v; 
    % Future changes
    load([''] % load data

    
    LEADTIME = []; 
    for i = 1:length(vstore100);
         xd = vv(i,1:300);
         xd(xd > 1) = 100; 
         kop = find(xd == 100) % & lat1002(i,1:300) > 24.5);
         if ~isempty(kop); % made landfall
             a = find(RIN2(i,1:kop(1)) > 45); % before landfall, RI
             if ~isempty(a); 
                 LEADTIME = [LEADTIME; kop(1) - a(1)];
             end
         end
    end

     cd('')
     save(['Rapid_intes_' SID{op} '_leadtime_25.mat'],'LEADTIME');   % save lead times from RI to landfall 
end

### Save data for plotting
close all
LEADTIMESSSP=[]
for i = 1:6;
    load(['Rapid_intes_' SID{i} '_leadtime.mat'])
    LEADTIME_A = [LEADTIME_A;LEADTIME];
    ncount = histc(LEADTIME,0:4:200);
    all_mean = ncount ./ nansum(ncount);
    LEADTIMESSSP(i,:) = all_mean;
end
SSP = (LEADTIMESSSP);

## Load for plotting 
load(['Rapid_intes_ncep_leadtime.mat']);
ncount = histc(LEADTIME,0:200);
all_mean = ncount ./ nansum(ncount);
nanmean(LEADTIME)

load(['Rapid_intes_ncep_leadtime_25.mat']);
ncount = histc(LEADTIME,0:200);
all_mean = ncount ./ nansum(ncount);
nanmean(LEADTIME)

%%
load(['Rapid_intes_ncep_leadtime_45.mat']);
ncount = histc(LEADTIME,0:200);
all_mean = ncount ./ nansum(ncount);
nanmean(LEADTIME)

X = 0:4:200;


