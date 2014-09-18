  list = {'T','h','c','nsquared','S','ks'};
% list = {'T','h','c','nsquared'};
% list = {'T','S','nsquared','h'}

start = round(length(tavg)/length(tavg)); % time index for start
tend = round(length(tavg));
n = length(list);
switch n
    case 1,
        figr = 1; figc = 1;
    case 2,
        figr = 2; figc = 1;
    case 3
        figr = 3; figc = 1;
    case 4
        figr = 4; figc = 1;
    case {5,6}
        figr = 3; figc = 2;
    case {7,8}
        figr = 4; figc = 2;
end

[ro,co] = size(Tavg);
% if co>1e4;    
    Myr = 1;
    tplot = tavg/(1e6*3.16e7);
%     t = (1:co)*tavg'/(1e6*3.16e7);
% else
%     tplot = tavg/(1000*3.16e7);
% end
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
for i=1:n
    switch list{i}
        case 'T'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,Tavg(:,start:tend));
            title('Temperature (K)','fontsize',14);
        case 'S'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,Savg(:,start:tend));
            title('Salinity (molal)','fontsize',14);
        case 'c'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,cavg(:,start:tend));
            title('Convective activity','fontsize',14);
        case 'nsquared'
            de = mgso4_dens(Savg(1:end-1,:),Tavg(1:end-1,:),pp(1:end-1,:));
            pde = mgso4_pden(Savg(2:end,:),Tavg(2:end,:),pp(2:end,:),pp(1:end-1,:));
            nsquared = g.*(de-pde)./(dz.*de);
            
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,nsquared(:,start:tend));
            title('Stratification (N^2, s^{-2})','fontsize',14);
        case 'q'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,10.^qavg,[-.4 .1]);
            title('Heat flux (W/m^2)','fontsize',14);
        case 'h'
            hthck = subplot(figr,figc,i,'align');
            plot(tplot(start:tend),havg(start:tend)/1e3);
            set(gca,'xlim',[tplot(start) tplot(tend)]);%,'ylim',[0 100]);
            ylabel('Thickness (km)','fontsize',14);
            title('Ice thickness (km)','fontsize',14);
        case 'kt'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,ktavg(:,start:tend));
            title('Mixing (KT, m^2/s)','fontsize',14);
        case 'ks'
            ax(i) = subplot(figr,figc,i,'align');
            imagesc(tplot(start:tend),z,ksavg(:,start:tend));
            title('Mixing (KS, m^2/s)','fontsize',14);
    end
%    set(gca,'XLim',[0 600])
    if (i > n - figc)
        if Myr
            xlabel('Time (My)','fontsize',14);
        else           
            xlabel('Time (ky)','fontsize',14);
        end
    end
    if ((mod(i,figc)==1 || figc == 1) && ~strcmp(list{i},'h'))
        ylabel('Depth (km)','fontsize',14);
    end
    if ~strcmp(list{i},'h')        
        hc = colorbar;
        if strcmp(list{i},'T')
            set(hc,'YDir','reverse');
        end
    end
end
linkaxes([ax hthck],'x');