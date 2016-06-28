% ==============================================================================
% Purpose: Make PV maps on different isopycnals using MIMOC climatology
%
% Dhruv Balwada
% 15/5/2016
% ==============================================================================

clear all
close all

% files naames of T and S files
fnames = dir('../MIMOC/data/MIMOC_Z_*');

% load T and S data
for mon = 1 : 12
    fname = ['../MIMOC/data/' fnames(mon).name];
    
    lon   = ncread(fname,'LONGITUDE');
    lat   = ncread(fname,'LATITUDE');
    pres  = ncread(fname,'PRESSURE');
    absSalt(:,:,:,mon)  = ncread(fname,'ABSOLUTE_SALINITY');
    consTemp(:,:,:,mon) = ncread(fname,'CONSERVATIVE_TEMPERATURE');
end

% convert to other variables (Potential temp, Practical Salinity and
% Neutral density

potTemp  = NaN*absSalt;
pracSalt = potTemp;
neutDens = potTemp;
potDens0 = potTemp;
[latGrid,lonGrid,presGrid] = meshgrid(lat,lon,pres);

for mon = 1 : 12
    potTemp(:,:,:,mon)  = gsw_pt_from_CT(squeeze(absSalt(:,:,:,mon)), squeeze(consTemp(:,:,:,mon)));
    pracSalt(:,:,:,mon) = gsw_SP_from_SA(squeeze(absSalt(:,:,:,mon)), presGrid, lonGrid, latGrid);
    potDens0(:,:,:,mon) = gsw_sigma0(squeeze(absSalt(:,:,:,mon)), squeeze(consTemp(:,:,:,mon)));
    neutDens(:,:,:,mon) = gamma_GP_from_SP_pt(pracSalt(:,:,:,mon), potTemp(:,:,:,mon), presGrid, lonGrid, latGrid) ;
    
    disp(mon)
end

% file names of ML files
fnames = dir('../MIMOC/data/MIMOC_ML_*');

% load T,S and depth for mixed layers
for mon = 1 : 12
    fname = ['../MIMOC/data/' fnames(mon).name];
    
    absSalt_ML(:,:,mon)  = ncread(fname,'ABSOLUTE_SALINITY_MIXED_LAYER');
    consTemp_ML(:,:,mon) = ncread(fname,'CONSERVATIVE_TEMPERATURE_MIXED_LAYER');
    depth_ML(:,:,mon)    = ncread(fname,'DEPTH_MIXED_LAYER');
end

% Estimate neutral density etc for mixed layer

[latGrid2d,lonGrid2d] = meshgrid(lat,lon);

for mon = 1 : 12
    potTemp_ML(:,:,mon)  = gsw_pt_from_CT(squeeze(absSalt_ML(:,:,mon)), squeeze(consTemp_ML(:,:,mon)));
    pracSalt_ML(:,:,mon) = gsw_SP_from_SA(squeeze(absSalt_ML(:,:,mon)), 0*lonGrid2d, lonGrid2d, latGrid2d);
    potDens0_ML(:,:,mon) = gsw_sigma0(squeeze(absSalt_ML(:,:,mon)), squeeze(consTemp_ML(:,:,mon)));
    neutDens_ML(:,:,mon) = gamma_GP_from_SP_pt(pracSalt_ML(:,:,mon), potTemp_ML(:,:,mon), 0*lonGrid2d, lonGrid2d, latGrid2d) ;
end

%
load('/tank/users/balwada/ARGO/Argo_RG_Climatology/ssh_SO_mean.mat');

%% Plots of surface conditions

for mon = 1:12;
    figure(1)
    clf
    set(gca,'fontsize',16)
    m_proj('Miller','lon',[0 360], 'lat',[-80 -30])
    hold all
    m_contourf(lonGrid2d, latGrid2d, depth_ML(:,:,mon),[0:20:300],'edgecolor','none')
    %     [c,h] = m_contour(lonGrid2d, latGrid2d, neutDens(:,:,1,mon),[26:0.2:28],'edgecolor','k')
    
    m_contour(ssh_mean.lon, ssh_mean.lat, ssh_mean.hmean, [-100 20],'edgecolor',[1 1 0.9],'linewidth',1.5)
    
    m_contour(lonGrid2d, latGrid2d, neutDens_ML(:,:,mon),[27.5 27.7 27.9],'edgecolor','r','linewidth',1.5)
    m_contour(lonGrid2d, latGrid2d, neutDens(:,:,21,mon),[27.5 27.7 27.9],'edgecolor','k','linewidth',1.5)
    
    %     clabel(c,h)
    m_coast('patch',[0 0 0])
    
    caxis([0 250])
    h=colorbar;
    ylabel(h,'Mixed Layer Depth (m)')
    m_grid
    title(datestr(mon*30,'mmm'))
    figName = ['./figures/surf_ML_27_7_' num2str(mon) '.png'];
    saveas(gcf,figName,'png')
end
%% Vertical sections
lonSel = 360-130;
idx = find(lon<lonSel,1,'last');

figure
hold all
[c,h]= contourf(lat, pres, squeeze(neutDens(idx,:,:,mon))',[22 26:0.2:28]);
clabel(c,h)
contour(lat, pres, squeeze(neutDens(idx,:,:,mon))',[27.6 27.6],'edgecolor','w','linewidth',2)
caxis([26 28])
set(gca,'Ydir','reverse')
axis([-80 0 0 2000])

%% Calcualte PV at all grid points
N2 = nan(length(lon), length(lat), length(pres)-1,12);
presMid = N2;
for mon =1:12
    
    for i=1:length(lon)
        for j=1:length(lat)
            [N2(i,j,:,mon),presMid(i,j,:,mon)]= gsw_Nsquared(squeeze(absSalt(i,j,:,mon)),...
                squeeze(consTemp(i,j,:,mon)),squeeze(presGrid(i,j,:)),squeeze(latGrid(i,j,:)));
        end
        
    end
    disp(mon)
end

rho0 = 1026;
dRhodZ = N2/9.81*rho0; % drho/dz for PV

f3D = (2*(2*pi/3600/24)*sind(latGrid(:,:,1:end-1)));
f2D = (2*(2*pi/3600/24)*sind(latGrid2d));
potVortZ = N2;

for mon =1:12
    potVortZ(:,:,:,mon) = f3D/rho0.*squeeze(dRhodZ(:,:,:,mon));
end

%% Interpolate and find the depth of a density layer

densSel = [27.4 27.5 27.6 27.7 27.8 27.9 28];
depDens = nan(length(lon), length(lat), length(densSel),12);

for mon = 1:12
    for i =1 : length(lon)
        for j =1:length(lat)
            
            for k =1:length(densSel)
                idk = find(neutDens(i,j,:,mon)<= densSel(k),1,'last');
                % if layer below the maximum depth of data set
                if idk == length(pres)
                    depDens(i,j,k,mon) = NaN;
                    continue;
                end
                % entire water column is denser (shoaled)
                if isempty(idk)
                    depDens(i,j,k,mon) = 0;
                    continue;
                end
                % over a land point
                if isempty(idk)  && isempty(find(~isnan(neutDens(i,j,:)),1))
                    depDens(i,j,k,mon) = NaN;
                    continue;
                end
                % for all other cases interpolate
                depDens(i,j,k,mon) = pres(idk) + (densSel(k) - neutDens(i,j,idk))* ...
                    (pres(idk+1) - pres(idk))/(neutDens(i,j,idk+1)-neutDens(i,j,idk));
            end
        end
    end
    disp(mon)
end

%% Depth plots
A = 27.7;
idDens = find(densSel == A);
for mon =1:12
    figure
    clf
    set(gca,'fontsize',16)
    m_proj('Miller','lon',[0 360], 'lat',[-80 -30])
    
    hold all
    m_contourf(lonGrid2d, latGrid2d, depDens(:,:,idDens,mon),[0:40:1000],'edgecolor','none')
    m_contour(lonGrid2d, latGrid2d, depDens(:,:,idDens,mon),[10 10],'edgecolor','r')
    
    m_coast('patch',[0 0 0])
    
    caxis([0 1000])
    h=colorbar;
    ylabel(h,['Depth (m) of ' num2str(A)] )
    m_grid
    title(datestr(mon*30,'mmm'))
    figName = ['./figures/surf_ML_27_7_' num2str(mon) '.png'];
    saveas(gcf,figName,'png')
end

%% Thickness plots
A = 27.7;
idDens = find(densSel == A);
for mon =1:12
    figure
    clf
    set(gca,'fontsize',16)
    m_proj('Miller','lon',[0 360], 'lat',[-80 -30])
    
    hold all
    m_contourf(lonGrid2d, latGrid2d, depDens(:,:,idDens+1,mon)-depDens(:,:,idDens-1,mon),[0:40:500],'edgecolor','none')
    m_contour(lonGrid2d, latGrid2d, depDens(:,:,idDens,mon),[10 10],'edgecolor','r')
    m_contour(ssh_mean.lon, ssh_mean.lat, ssh_mean.hmean, [-100 20],'edgecolor','k','linewidth',1.5)
    
    %     clabel(c,h)
    m_coast('patch',[0 0 0])
    
    caxis([0 500])
    h=colorbar;
    ylabel(h,['Thickness (m) of ' num2str(A)] )
    m_grid
    title(datestr(mon*30,'mmm'))
    figName = ['./figures/thick_27_7_' num2str(mon) '.png'];
    saveas(gcf,figName,'png')
end

%% PV on layer plots
A = 27.7;
dDens = mean(diff(densSel));
idDens = find(densSel == 27.7);

for mon =1:12
    figure(1)
    clf
    set(gca,'fontsize',16)
    m_proj('Miller','lon',[0 360], 'lat',[-80 -30])
    
    pvPlot = f2D.*dDens/rho0./(depDens(:,:,idDens+1,mon)-depDens(:,:,idDens-1,mon));
    
    hold all
    m_contourf(lonGrid2d, latGrid2d, pvPlot,...
        [-10:0.1:0]*10^(-10),'edgecolor','none')
    m_contour(lonGrid2d, latGrid2d, pvPlot,...
        [-1:0.2:0]*10^(-10),'edgecolor','k')
    %     m_contour(lonGrid2d, latGrid2d, depDens(:,:,idDens,mon),[10 10],'edgecolor','r','linewidth',2)
    %     m_contour(lonGrid2d, latGrid2d, depDens(:,:,idDens,mon),[100 100],'edgecolor','b','linewidth',2)
    
    m_contour(lonGrid2d, latGrid2d, neutDens_ML(:,:,mon),[A-dDens A A+dDens],'edgecolor','r','linewidth',1.5)
    m_contour(lonGrid2d, latGrid2d, neutDens(:,:,21,mon),[A-dDens A A+dDens],'edgecolor','k','linewidth',1.5)
    m_contour(ssh_mean.lon, ssh_mean.lat, ssh_mean.hmean, [-100 20],'edgecolor',[0.5 0.5 0.9],'linewidth',1.5)
    
    %     clabel(c,h)
    m_coast('patch',[0 0 0])
    
    caxis([-8 0]*10^(-10))
    h=colorbar;
    ylabel(h,['Thickness (m) of ' num2str(A)] )
    m_grid
    title(datestr(mon*30,'mmm'))
    figName = ['./figures/thick_27_7_' num2str(mon) '.png'];
    saveas(gcf,figName,'png')
end
