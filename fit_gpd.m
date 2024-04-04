% function for selecting threshold and fitting gpd distribution 

function [fit, u, type] = fit_gpd(data, minp)
data = data(data~=Inf);
data = sort(data);
n = length(data);
myfun1 = @(x)ksdensity(x,'Bandwidth',4,'Function','cdf');
%if 0.2*n<minp
    fit = fitdist(data', 'Kernel', 'Bandwidth', 4); 
    u = 0; 
    type = 'kde'; 
%
end

