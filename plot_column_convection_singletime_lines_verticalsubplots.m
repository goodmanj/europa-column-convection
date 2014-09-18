  list = {'T','S','c','ks','nsquared'};%,'Rrho'};
% list = {'T','h','c','nsquared'};
% list = {'T','S','nsquared','h'}
t_inds = [1 39];

start = round(length(tavg)/length(tavg)); % time index for start
tend = round(length(tavg));
n = length(list);
figc = 1;
figr = n;

Myr = 0; kyr = 0;
[ro,co] = size(Tavg);
 if tavg(end)>1e6*3.16e7;    
     Myr =1;
    tplot = tavg/(1e6*3.16e7);
     t = (1:co)*tavg'/(1e6*3.16e7);
 elseif tavg(end)>10e3*3.16e7
     kyr = 1;
     tplot = tavg/(1000*3.16e7);
 else
     tplot = tavg/3.16e7;
 end
    z = ((ro:-1:1)-.5)*dz/1000;  

[hgrid,zgrid] = meshgrid(havg,((Nz-1):-1:0)*dz);
pp = ((rhoi*hgrid*g) + rhow*g*zgrid)*patobar;

%switch the x and y orientation of the plot matrix
% holdr = figr ;
% figr = figc;
% figc = holdr;

S0str = num2str(S0);S0str(S0str=='.')='p';
Sfstr = num2str(Sfluxseafloor); Sfstr(Sfstr=='-')='m';
Hstr = num2str(Hseafloor); Hstr(Hstr=='.')='p';
Tstr = num2str(t/3.16e13); Tstr(Tstr=='.') = 'p';
hf = figure;set(hf,'name',['S0_' S0str 'molal_T' Tstr 'Myr_H' Hstr 'Wm2_Sf' Sfstr ]);
for is=1:n
    switch list{is}
        case 'T'
            ax(is) = subplot(figr,figc,is,'align');
            hold all
            h_pl = plot(z,Tavg(:,t_inds));
            ylabel('Temperature (K)','fontsize',14);
            box on
        case 'S'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z,Savg(:,t_inds));
            ylabel('Salinity (molal)','fontsize',14);
        case 'c'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z(2:end),cavg(:,t_inds));
            ylabel('Convective activity','fontsize',14);
        case 'nsquared'
            de = mgso4_dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            pde = mgso4_pden(Savg(2:end,:),Tavg(2:end,:),pp(2:end,:),pp(1:end-1,:));
            nsquared = g.*(de-pde)./(dz.*de);
            
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z(2:end),nsquared(:,t_inds));
            ylabel('N^2 (s^{-2})','fontsize',14);
        case 'q'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z,10.^qavg(:,t_inds));
            ylabel('Heat flux (W/m^2)','fontsize',14);
        case 'h'
            hthck = subplot(figr,figc,is,'align');
            h_pl = plot(z,havg(t_inds)/1e3);
            set(gca,'xlim',[tplot(start) tplot(tend)]);%,'ylim',[0 100]);
            ylabel('Ice thickness (km)','fontsize',14);
        case 'kt'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z,ktavg(:,t_inds));
            ylabel('K_T (m^2/s)','fontsize',14);
        case 'ks'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z(2:end),ksavg(:,t_inds));
            ylabel('K_S (m^2/s)','fontsize',14);
        case 'Rrho'
            ax(is) = subplot(figr,figc,is,'align');
            h_pl = plot(z(2:end),Rrho(:,t_inds));
            ylabel('R_{\rho}','fontsize',14);
    end
end
xlabel('z_{ocean} (km)','fontsize',14);

if Myr
    tstr = 'Time (My)';
elseif kyr           
    tstr = 'Time (ky)';
else
    tstr = 'Time (yr)';           
end

l_str = cell(1,length(t_inds));
for it = 1:length(t_inds)
    l_str{it} = [tstr ' : ' num2str(tplot(t_inds(it)))];
end
    legend(h_pl,l_str)
linkaxes(ax,'x');