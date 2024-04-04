%% RI calculation example
$ length is number of events, example

# Load events
load('/Users/joelockwood/PhD/Research/Covarianceproject/ADCIRC_cmip6/StormFiles/histo_storms/UScoast6_AL_ecearth6_20thcal_roEst1rmEst1_trk100.mat')
for i = 1:length ;
    vmax = vstore100(i,:);
    vstore1002(i,1:length(xi)) = lat10;
end

% Rapid intensificition, vmax over time
 
for i=1:length;
    asd = ~isnan(vstore1002(i,:));
    for id=1:length(asd)-time;
        RIN(i,id) = vstore1002(i,id+time) - vstore1002(i,id); 
    end
    RINlt(i) = nanmax(RIN(i,:));
end


%%%%%% Percent undergoing RI %%%%%%%%

RIN2n =[]; RINn =[]; 
for v = 1:size;
    if any(RIN2(v,:) > 35); RIN2n(v) = 1; end
    if any(RIN(v,:) >  35); RINn(v) = 1; end 
end

%%%%%%%%%% Plotting maps %%%%%%%
m_proj('equidist','lat',[20 45],'lon',[-100 -60],'sphere','aspect',.8);
set(gca,'linewidth',.2)
close all
lon1002(lon1002 == 0) = nan;
lat1002(lat1002 == 0) = nan;

for v = 1:time;
    if RINn(v) > 0 ; hold on; adss = find(RIN2(v,1:end) > 35):find(RIN2(v,1:end)> 35); m_scatter(lon1002(v,:), lat1002(v,:), [],RIN2(v,adss),'filled')   ; end
end
m_coast('color',[.5 .5 .5],'linewi',2);

 
%%%%%%%% Rapid intensification maps example  %%%%%%%%%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.04], [0.1 0.07], [0.05 0.1]);
if ~make_it_tight,  clear subplot;  end

ctype = 'div'
cname = 'RdBu'
ncol = 100
[colormap3]=cbrewer(ctype, cname, ncol)

close all;
m_proj('equidist','lat',[0 45],'lon',[-100 -60],'sphere','aspect',.8);
hold on;

subplot(1,3,1)
m_pcolor(u,v,nanmean(meantrans1,3).')  % plot a map of the RIR
m_grid('linest','none','xtick',16,'tickdir','in','xtick',[],'ytick',[]);
ctype = 'div'
cname = 'RdBu'
ncol = 100;
[colormap3]=cbrewer(ctype, cname, ncol); hold on; 
colormap((jet)); colorbar()
hold on;
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
    m_plot(y,vv,'color',[.5 .5 .5],'linewidth',2)
end

