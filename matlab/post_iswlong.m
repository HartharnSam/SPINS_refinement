% Post processing and plotting for u_top and u_bottom.

clearvars

params = spins_params();
Nx = params.Nx;
t_slice = params.plot_interval/params.plot_interval_1d;
t_full  = params.final_time   /params.plot_interval;
t_all   = params.final_time   /params.plot_interval_1d;

u_top = zeros(Nx,t_all);
u_bottom = zeros(Nx,t_all);

for ii = 1:t_full
    u_top(:,(1:t_slice)+(ii-1)*t_slice)=u_top_reader(ii);
    u_bottom(:,(1:t_slice)+(ii-1)*t_slice)=u_bottom_reader(ii);
end

x2d = xgrid_reader;
x1d = x2d(:,1);
t1d = linspace(0,params.final_time,t_all);
[xx,tt] = meshgrid(x1d,t1d);

set(0,'defaulttextinterpreter','latex')
f1 = figure('visible','off'); clf, colormap temperature
set(gcf,'units','inches')
set(gcf,'position',[1 1 7 7])
pcolor(xx,tt,u_top'), shading flat
caxis(max(abs(u_top(:)))*[-1,1])
title('u\_top'), xlabel('$x$ (m)'), ylabel('$t$ (s)')
saveas(f1,'u_top','png')

set(0,'defaulttextinterpreter','latex')
f2 = figure('visible','off'); clf, colormap temperature
set(gcf,'units','inches')
set(gcf,'position',[1 1 7 7])
pcolor(xx,tt,u_bottom'), shading flat
caxis(max(abs(u_bottom(:)))*[-1,1])
title('u\_bottom'), xlabel('$x$ (m)'), ylabel('$t$ (s)')
saveas(f2,'u_bottom','png')

quit
