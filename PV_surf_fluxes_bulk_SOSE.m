% ==============================================================================
% Purpose : Estimates the surface PV fluxes
% Uses bulk formulae listed in Deremble et al 2014
% Dhruv Balwada
% 20/6/2016
% ==============================================================================
clear all 
close all 

dirData  = '/tank/groups/SO/SOSE_data/';
fileFW   = 'oceFWflx.0000000100'; 
fileQnet = 'TFLUX.0000000100'; 
fileSalt = 'SALT.0000000100'; 
fileTemp = 'THETA.0000000100'; 
fileGrid = 'grid.mat';
fileDens = 'GAMMA.0000000100'; 
fileTaux = 'oceTAUX.0000000100'; 
fileTauy = 'oceTAUY.0000000100'; 
fileMLD  = 'MLD.0000000100';
gridSOSE = load([dirData fileGrid]);

minrec = 1;
maxrec = 438; 

FWflux = nan(size(gridSOSE.XC,1), size(gridSOSE.XC,2), maxrec-minrec+1);
Qnet   = FWflux;
tempS  = nan(size(gridSOSE.XC,1), size(gridSOSE.XC,2), length(gridSOSE.RC));
tempT  = tempS;
tempD  = tempS;
SSS    = FWflux; 
SST    = FWflux; 
SSD    = FWflux; % surface density
tauX   = FWflux;
tauY   = FWflux;
MLD    = nan(size(gridSOSE.XC,1), size(gridSOSE.XC,2), 5*maxrec-minrec+2); 
MLD_5dmean = FWflux; 

for n = minrec:maxrec
    FWflux(:,:,n) = rdmds([dirData fileFW],NaN,'rec',n);
    Qnet(:,:,n)   = rdmds([dirData fileQnet],NaN,'rec',n);
    tauX(:,:,n)   = rdmds([dirData fileTaux],NaN,'rec',n);
    tauY(:,:,n)   = rdmds([dirData fileTauy],NaN,'rec',n);
    tempS         = rdmds([dirData fileSalt],NaN,'rec',n);
    tempT         = rdmds([dirData fileTemp],NaN,'rec',n);
    tempD         = rdmds([dirData fileDens],NaN,'rec',n);
    SSS(:,:,n)    = tempS(:,:,1);
    SST(:,:,n)    = tempT(:,:,1);
    SSD(:,:,n)    = tempD(:,:,1);
    
    disp(n)
end

% MLD file currently is a daily output 
for n = minrec:maxrec*5+1
    MLD(:,:,n)    = rdmds([dirData fileMLD],NaN,'rec',n);
end

% Average daily outputs to 5 day.
p=1;
for n = minrec:maxrec
    MLD_5dmean(:,:,n) = nanmean(MLD(:,:,1+(n-1)*5:n*5),3);
end

% Parameters 
alpha = gsw_alpha(SSS, SST, 0); 
beta  = gsw_beta(SSS, SST, 0);  

cp    = gsw_cp_t_exact(SSS, SST, 0);
f     = 2*2*pi/24/3600*sind(gridSOSE.YC); 

rho0  = 1026;

% Diabatic J 
JzBQnet = nan*FWflux; 
JzBEP   = JzBQnet; 

for n = minrec:maxrec
    JzBQnet(:,:,n) = -f./MLD_5dmean(:,:,n).*alpha(:,:,n).*Qnet(:,:,n)./cp(:,:,n); 
    % There is minus sign here as E-P = - (FWflux) 
    JzBEP(:,:,n)   = f.*beta(:,:,n).*SSS(:,:,n).*FWflux(:,:,n)./MLD_5dmean(:,:,n);
    disp(n)
end

% Mechanical J 

JzF = nan*FWflux; 
dOmdx = nan*FWflux(:,:,1); 
dOmdy = nan*FWflux(:,:,1); 

for n = minrec:maxrec
    for nx = 2:size(gridSOSE.XC,1)-1
        for ny = 2:size(gridSOSE.XC,2)-1
            dOmdx(nx,ny) = (SSD(nx+1, ny,n) - SSD(nx-1, ny,n))./gridSOSE.DXC(nx,ny)/2;
            dOmdy(nx,ny) = (SSD(nx, ny+1,n) - SSD(nx, ny-1,n))./gridSOSE.DYC(nx,ny)/2;
        end
    end
    % add the part about derivatives at the periodic boundary
    
    disp(n)

    JzF(:,:,n) = (tauX(:,:,n).*dOmdy - tauY(:,:,n).*dOmdx) ...
                /rho0/0.4.*f./sqrt(sqrt(tauX(:,:,n).^2 + tauY(:,:,n).^2)/rho0);
%    disp(n)
end

% J from entrainement 

JzW = nan*FWflux; 
g =9.81; 
for n=minrec:maxrec
    JzW(:,:,n) =  0.7/g./(MLD_5dmean(:,:,n).^2).*f*rho0.*(sqrt(tauX(:,:,n).^2 + tauY(:,:,n).^2)/rho0).^(3/2);
end

% Calculate Means

% Eulerian means (no conditional averaging on outcrops)
% provides a sense of PV flux into the ocean
% these are the kind of calculations that were done in Czaja & Haussmann 2009
emeanJzF     = nanmean(JzF,3); 
emeanJzW     = nanmean(JzW,3); 
emeanJzBEP   = nanmean(JzBEP,3);
emeanJzBQnet = nanmean(JzBQnet,3);
emeanSSD     = nanmean(SSD,3);


% Plots
plot_flag = 0;

if plot_flag == 1
num = 7;
close all 
figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzF,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z^F')
saveas(gcf,'JZF_emean.png','png')
%

figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzW,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z^W')
saveas(gcf,'JZW_emean.png','png')

%

figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzBEP,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z^B - EP')
saveas(gcf,'JZB_EP_emean.png','png')


%

figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzBQnet,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z^B - Qnet')
saveas(gcf,'JZB_Qnet_emean.png','png')

%

figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzBQnet + emeanJzBEP,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z^B - Qnet + EP')

saveas(gcf,'JZB_emean.png','png')

%

figure 
colormap(redbluecmap)
set(gca,'fontsize',16)
m_proj('Miller','lon',[0 360], 'lat', [-80 -25])
m_contourf(gridSOSE.XC, gridSOSE.YC, emeanJzBQnet + emeanJzBEP + emeanJzF + emeanJzW,[-1000 -num:num]*10^-12,'edgecolor','none')
hold all
% m_contour(gridSOSE.XC, gridSOSE.YC, emeanJzF,[0 0],'edgecolor','k')
[c,h] = m_contour(gridSOSE.XC, gridSOSE.YC, emeanSSD, [25:0.2:28],'edgecolor','k');
clabel(c,h)
caxis([-num num]*10^(-12))
m_coast('patch',[ 0 0 0])
m_grid
colorbar
title('J_z')

saveas(gcf,'JZ_emean.png','png')

end
