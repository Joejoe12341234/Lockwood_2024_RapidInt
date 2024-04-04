clear all;
close all; m_proj('equidist','lat',[20 50],'lon',[-99 -60],'sphere','aspect',10);
cd('/Users/joelockwood/PhD/Research/Rapid_intense_project/data')
file = 'ETOPO_2022_v1_60s_N90W180_bed.nc'
ncdisp(file)
lat = ncread(file, 'lat');
lon = ncread(file, 'lon');
f=find(lat>18);; f=f(1);
f2=find(lat<55);;f2=f2(end);
l=find(lon<-100);;l=l(end);
l2=find(lon>-70);;l2=l2(1);
z = ncread(file, 'z');

f=find(lat>18);; f=f(1);
f2=find(lat<55);;f2=f2(end);
l=find(lon<-100);;l=l(end);
l2=find(lon>-70);;l2=l2(1);
lat2 = lat(f:15:f2);
lon2 = lon(l:15:l2);
z2=z(l:15:l2,f:15:f2).';
SID={'canesm','cnrm6','ecearth6','ipsl6','miroc6','ukmo6'}

for v = 1; 
    op = v; 
    % Future changes
    load(['/Users/joelockwood/PhD/Research/Rapid_intense_project/data/ncep_storms/UScoast6_AL_ncep_reanal_roEst1rmEst1_trk100.mat'])
    % load(['/Users/joelockwood/PhD/Research/Covarianceproject/ADCIRC_cmip6/StormFiles/ssp245/UScoast6_AL_'   SID{v}  '_ssp245cal_roEst1rmEst1_trk100.mat'])

    lon1002 =[]; lat1002 =[] 
    for i = 1:length(vstore100) ;
        vmax = vstore100(i,:);
        vmax(vmax == 0) = nan;
        asd = ~isnan(vmax);
        times2 = size(vmax(asd),2);
        time = 1:times2;
        xi = 1:0.5:times2;
        o = vstore100(i,asd);   
        lat10 = interp1(time,o,xi,'linear');
        vstore1002ft(i,1:length(xi)) = lat10;
        o = lat100(i,asd);
        lat10 = interp1(time,o,xi,'linear');
        lat1002(i,1:length(xi)) = lat10;
        o = lon100(i,asd);
        lat10 = interp1(time,o,xi,'linear');
        lon1002(i,1:length(xi)) = lat10;
    end

    vv = zeros(size(vstore1002ft,1),size(vstore1002ft,2)); % dummy bathymetry
    for i = 1:length(vstore100) 
        for o = 1:size(lon1002,2);
            [ d, loni ] = min(abs(lon1002(i,o)-lon2 ) );
            [ d, lati ] = min( abs(lat1002(i,o)-lat2 ) );
            if z2(lati,loni) > 0; 
                vv(i,o) = z2(lati,loni);
            else
                vv(i,o) = z2(lati,loni);
            end;
            if vstore1002ft(i,o) == 0 ; vv(i,o) = nan; end
        end
    end

    lgnt = 10;
%  %   Find differences between adjacent elements
%     for i = 1:length(vstore100);
%        xd = vv(i,:);
%         xd(xd > 1) = 100; 
%         vstore1002ft2 = vstore1002ft;
%         % subplot(2,1,1); hold on; m_scatter(lon1002(i,1:400),lat1002(i,1:400),3,'k','filled')
% %         for v = 40:350;
% %             if xd(v-5) == 100 & lat1002(i,v) > 36;
% %                vstore1002ft(i,v-5:end) = nan;
% %                lat1002(i,v-5:end) = nan;
% %                lon1002(i,v-5:end) = nan;
% %                % disp(v)
% %                break;
% %             end
% %         end
% % 	    g = find(xd == 100); 
% %         vstore1002ft(i,g) = nan;
% %         lat1002(i,g) = nan;
% %         lon1002(i,g) = nan;
%    end
    
    vstore1002ft(vstore1002ft == 0) = nan;
    vstore1002grad=[]; RIN2=zeros(size(lon1002,1),size(lon1002,2));
    for i=1:length(vstore100); % calculate the RI 
         asd = ~isnan(vstore1002ft(i,:));
         for id=1:size(lon1002,2)-24;
              RIN2(i,id) = vstore1002ft(i,id+23) - vstore1002ft(i,id); 
         end
         RIN2lt(i) = nanmax(RIN2(i,:));
    end
    
    lon1002f = lon1002;
    lat1002f = lat1002;
    RIN2n =[]; RINn =[];
    
    for vv = 1:length(vstore100); % if any have RI above 30
        if any(RIN2(vv,:) > 45); RIN2n(vv) = 1; end
    end
 
      cd('/Users/joelockwood/PhD/Research/Rapid_intense_project/Saving/new/')
     %save(['Rapid_intes_ncep.mat'],'RIN2n','RIN2lt','RIN2','lon1002f','lat1002f','vstore1002ft');    
     save(['Rapid_intes_' SID{op} '_45.mat'],'RIN2n','RIN2lt','RIN2','lon1002f','lat1002f','vstore1002ft');    

     latTC = round(lat1002f); lonTC = round(lon1002f);rr = RIN2 ; rr(isnan(rr)) = [];
     rr = reshape(rr,[],1) ;asd = find(rr > 45); asd2 = find(rr > -100); latTC = reshape(latTC ,[],1) ;
     lonTC = reshape(lonTC ,[],1) ;latTCri = latTC(asd); lonTCri = lonTC(asd); latTC = latTC(asd2); lonTC = lonTC(asd2);
     ssp_all=[]; ssp_ri=[];
            for u = -120:-0;
                for v = 1:55;
                     aa = (find(latTCri == v & lonTCri == u));
                     ssp_ri(u+121,v) = length(aa);
                     aa = (find(latTC == v & lonTC == u));
                     ssp_all(u+121,v) = length(aa);
                end
            end
      cd('/Users/joelockwood/PhD/Research/Rapid_intense_project/Saving/new/')
      ratioo = [];ssp_RIR=[];
      ssp_RIR = (ssp_ri./(ssp_all)).';
      ratioo = (ssp_ri./(ssp_all)).';
     % save(['Rapid_intes_' SID{op} '_ratio_25.mat'],'ssp_RIR');    
      save(['Rapid_intes_ncep_ratio_45.mat'],'ratioo');    
     % save(['Rapid_intes_' SID{op} '_ratio.mat'],'ratioo','ncep_ri','ncep_all');   
end;

temptrag2=[];
for q = 1:length(SID);
    FileName   = 'temp_diff.mat';
    D = ['/Users/joelockwood/PhD/Research/Rapid_intense_project/data/temp/' SID{q} '_temp_diff.mat' ];
    File  = fullfile(D, FileName);
    load(D); 
    tempdiffmap = nanmean(nanmean(variable_mean,2),1);
    temptrag2(q) = tempdiffmap;
end

%%
figure(2)
cd('/Users/joelockwood/PhD/Research/Rapid_intense_project/Saving/new/')
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.0001 0.025], [0.07 0.07], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end
u = -120:-0
v = 1:55
set(0, 'DefaultAxesFontSize', 14);
set(0,'DefaultAxesTitleFontWeight','normal');
set(0, 'defaultAxesFontName', 'Arial');
close 
close all;

load('Rapid_intes_ncep_ratio.mat')
f = subplot(1,3,1)
m_proj('equidist','lat',[23 44.2],'lon',[-99 -72],'sphere','aspect',10);
hold on;
% savingVMAXncep_ratio(savingVMAXncep_ratio == 0) = nan;
ratioo(ratioo == 0) = nan;
A1 = ratioo;
A1 = smooth2a(A1,2,2);
m_contourf(u,v,A1,35,'edgecolor','none') 
caxis([0 0.1])
text( 0.0,1.06,'a) NCEP rapid intensification ratio (RIR)','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[10 15 20 25 30 35 40 45]);
ctype = 'seq'
cname = 'OrRd'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap((parula)); c=colorbar('southoutside');
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end
text( 0.0,1.06,'a) NCEP rapid intensification ratio (RIR)','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[10 15 20 25 30 35 40 45]);
ctype = 'seq'
cname = 'OrRd'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(f,(colormap3)); c=colorbar('southoutside'); caxis([-0.1 0.1])
caxis([0.02 0.07])

file = '/Users/joelockwood/PhD/Research/Rapid_intense_project/IBTrACS.ALL.v04r00-2.nc'
 ncdisp(file)
usa_wind = ncread(file,'usa_wind');
lat = ncread(file,'lat').';
lon = ncread(file,'lon').'; name = ncread(file,'name'); usa_wind = usa_wind.';
time  = ncread(file,'time'); 
 basin = ncread(file,'basin'); dt = datetime(time.*24.*3600, 'ConvertFrom', 'epochtime', 'Epoch', '1858-11-17');
 YEARS = dt.Year.' ;MONTHS = dt.Month;DAYS = dt.Day;HOURS = dt.Hour;
for v = 1:13621;
  if any(YEARS(v,:) > 1990) ; 
       if any(YEARS(v,:) < 2010);
            kkko(v) = 1;
        else 
            kkko(v) = 0;
        end
       else 
      kkko(v) = 0;
     end
 end
% 
sop = find((kkko == 0));
 lat(sop,:) =[];
lon(sop,:) =[];
 usa_wind(sop,:) =[];
 RIN2=[];RIN2lt=[];
for i = 1:length(usa_wind);
     asd = ~isnan(usa_wind(i,:));
     for id=1:size(usa_wind,2)-9;
        RIN2(i,id) = (usa_wind(i,id+8) - usa_wind(i,id)); 
    end
    RIN2lt(i) = nanmax(RIN2(i,:));
end
RIN2h=RIN2;

% % 3 hourly timesteps 
K=1; namez=[];
for v = 1:length(RIN2lt)
    if any(RIN2h(v,:) > 22);
        a = find(RIN2h(v,:) > 22);
        hold on; m_scatter(lon(v,a(1,1)),lat(v,a(1,1)),45,'k','filled')
         namez(:,K) = name(:,v);
         K=K+1
     end
 end
 hold on; K=1; 
 for v = 1:length(RIN2lt)
     if any(RIN2h(v,:) > 30);
         a = find(RIN2h(v,:) > 30);
         for vv = 1:length(a);
         hold on; m_scatter(lon(v,a(1,1)),lat(v,a(1,1)),45,'k','filled')
        end
        namez(:,K) = name(:,v);
         K=K+1
     end;
 end;
 
SSO = [];
for i = 1:6; 
    load(['Rapid_intes_' SID{i} '_ratio.mat']);    
    SSO(:,:,i) = ssp_RIR;
end;

u = -120:-0;
v = 1:55;
subplot(1,3,2);
u = -120:-0
v = 1:55
A = 100.*((nanmean(SSO(:,:,:),3) - ratioo) ./ ratioo);
A = smooth2a(A,2,2);
m_contourf(u,v,A,35,'edgecolor','none') 
text( 0.0,1.06,'b) Percent change in RIR','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[]);
ctype = 'div'
cname = 'RdBu'
ncol = 100;
c=colorbar('southoutside'); 
ylabel(c,'%','FontSize',15); 
ctype = 'div'
cname = 'RdBu'
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(flipud(colormap3)); 
caxis([40 150])
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end

subplot(1,3,3)
u = -120:-0
v = 1:55
hold on;
for vvv = 1:6;
    AA(:,:,vvv) = (100.*((SSO(:,:,vvv) - ratioo) ./ ratioo)) ./ temptrag2(vvv);
end
A = smooth2a(nanmean(AA,3),2,2);
m_contourf(u,v,A,35,'edgecolor','none') 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[]);
hold on; 
text( 0.0,1.06,'c) Percent change in RIR per degree global warming','Units','Normalized', 'Color','k','FontSize',20); 
ctype = 'div'
cname = 'RdBu'
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(flipud(colormap3)); 
c=colorbar('southoutside'); 
caxis([0 50])
ylabel(c,'% per global degree warming','FontSize',15); 
w = c.LineWidth;
c.LineWidth = 1.2;
hold on;
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end

ctype = 'seq'
cname = 'YlOrRd'
ncol = 200
[colormap2]=cbrewer(ctype, cname, ncol);
colormap(f,colormap2)

%%

figure(2)
cd('/Users/joelockwood/PhD/Research/Rapid_intense_project/Saving/new/')
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.0001 0.025], [0.07 0.07], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end
u = -120:-0
v = 1:55
set(0, 'DefaultAxesFontSize', 14);
set(0,'DefaultAxesTitleFontWeight','normal');
set(0, 'defaultAxesFontName', 'Arial');
close 
close all;

load('Rapid_intes_ncep_ratio_45.mat')
f = subplot(1,3,1)
m_proj('equidist','lat',[23 44.2],'lon',[-99 -72],'sphere','aspect',10);
hold on;
ratioo(ratioo == 0) = nan;
A1 = ratioo.*2;
A1 = smooth2a(A1,2,2);
m_contourf(u,v,A1,35,'edgecolor','none') 
caxis([0 0.1])
text( 0.0,1.06,'a) NCEP rapid intensification ratio (RIR)','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[10 15 20 25 30 35 40 45]);
ctype = 'seq'
cname = 'OrRd'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap((parula)); c=colorbar('southoutside');
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end
text( 0.0,1.06,'a) NCEP rapid intensification ratio (RIR)','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[10 15 20 25 30 35 40 45]);
ctype = 'seq'
cname = 'OrRd'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(f,(colormap3)); c=colorbar('southoutside'); caxis([-0.1 0.1])
caxis([0.0 0.04])

SSO = [];
for i = 1:6; 
    load(['Rapid_intes_' SID{i} '_ratio_45.mat']);    
    SSO(:,:,i) = ssp_RIR;
end;
u = -120:-0;
v = 1:55;
subplot(1,3,2);
u = -120:-0
v = 1:55
A = 100.*((nanmean(SSO(:,:,:),3) - ratioo) ./ ratioo);
A = smooth2a(A,2,2);
m_contourf(u,v,A,35,'edgecolor','none') 
text( 0.0,1.06,'b) Percent change in RIR','Units','Normalized', 'Color','k','FontSize',20); 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[]);
ctype = 'div'
cname = 'RdBu'
ncol = 100;
c=colorbar('southoutside'); 
ylabel(c,'%','FontSize',15); 
ctype = 'div'
cname = 'RdBu'
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(flipud(colormap3)); 
caxis([40 250])
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end

subplot(1,3,3)
u = -120:-0
v = 1:55
hold on;
for vvv = 1:6;
    AA(:,:,vvv) = (100.*((SSO(:,:,vvv) - ratioo) ./ ratioo)) ./ temptrag2(vvv);
end
A = smooth2a(nanmean(AA,3),2,2);
m_contourf(u,v,A,35,'edgecolor','none') 
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[-90 -80 -70 -60 -30],'ytick',[]);
hold on; 
text(0.0,1.06,'c) Percent change in RIR per degree global warming','Units','Normalized', 'Color','k','FontSize',20); 
ctype = 'div'
cname = 'RdBu'
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap(flipud(colormap3)); 
c=colorbar('southoutside'); 
caxis([0 150])
ylabel(c,'% per global degree warming','FontSize',15); 
w = c.LineWidth;
c.LineWidth = 1.2;
hold on;
load('/Users/joelockwood/PhD/Research/Rapid_intense_project/MIROC6_future/USstates.mat')
k = 1:52;
k(26) = [];
k(48) = [];
k(8) = [];
m_coast('color',[0.3 .3 .3],'linewi',3); hold on
m_coast('patch',[.5 .5 .5],'edgecolor','none');
for rr = 1:49
    rr = k(rr)
    vv = state_y{rr}
    y = state_x{rr}  
    hold on; 
    m_plot(y,vv,'color',[0.3 .3 .3],'linewidth',2)
end

ctype = 'seq'
cname = 'OrRd'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap((colormap3));

ctype = 'seq'
cname = 'YlOrRd'
ncol = 200
[colormap2]=cbrewer(ctype, cname, ncol);
colormap(f,colormap2)
