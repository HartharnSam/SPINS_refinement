%DRIVER_ISWLONG     u = w = 0, rho = rhobar(z - eta), where eta has both a
% seiche and a sech^2 disturbance. No interpolation needed.

clearvars

test = false;                       % set to false to write data to disk

% Case independent parameters
f  = 0.00;
g  = 9.81;
nu = 1e-6;
kappa = 2e-7;
pert  = 0.00;
final_time = 150;
plot_interval = 1;
plot_interval_1d = 0.0;
slip = true;

% Case dependent parameters
par = readtable('tilttank.txt');     % this may not work for R2017b or later

%% Loop through different cases
for numcase = 1:height(par)
    casename = char(par.casename(numcase));

    % Grid
    Lx = par.Lx(numcase);       % this determines the width of the seiche
    Lz = par.Lz(numcase);
    Nx = par.Nx(numcase);
    Nz = par.Nz(numcase);
        
    x = Lx*(0.5:Nx-0.5)/Nx;     % free-slip in x
%    z = Lz*(0.5:Nz-0.5)/Nz;     % free-slip in z
    [~,ztmp] = cheb(Nz);            % no-slip in z: size(z1d)==Nz+1
    z = flipud(0.5*Lz*(ztmp+1));    % so that z is in increasing order
    [xx,zz] = meshgrid(x,z);	% 2D grid in Matlab's meshgrid format
    
    % Stratification
    drho = 0.01;                % density difference
    z0   = par.z0(numcase);     % pycnocline location
    d    = par.d(numcase);      % pycnocline thickness
    amp  = par.amp(numcase);	% amplitude of mode-1 wave
    wl   = par.wl(numcase);     % width of mode-1 wave
    tilting = par.tilting(numcase);	% slope of the tank

    rho = 1 - 0.5*drho*tanh((zz - z0 + xx*tilting + ...
        amp*exp(-(xx/wl).^2))/d);
    
    % Velocity
    u = zeros(size(rho));
    w = zeros(size(rho));

    if test     % plot
        figure(numcase)
        pcolor(xx,zz,rho), shading flat, colormap darkjet
        
    else        % write data to disk
        mkdir(['../' casename]), cd(['../' casename])
        
        % SPINS grid ordering
        dspins = permute(rho,[2,1]);
        uspins = permute(u,[2,1]);
        wspins = permute(w,[2,1]);
        
        fid = fopen('rho.orig','wb'); fwrite(fid,dspins,'double'); fclose(fid);
        fid = fopen('u.orig','wb'); fwrite(fid,uspins,'double'); fclose(fid);
        fid = fopen('w.orig','wb'); fwrite(fid,wspins,'double'); fclose(fid);
        
        % Write parameters to spins.conf
        fid = fopen('spins.conf','wt');
        fprintf(fid,'Lx = %12.8f \n',Lx);
        fprintf(fid,'Lz = %12.8f \n',Lz);
        fprintf(fid,'Nx = %d \n',Nx);
        fprintf(fid,'Nz = %d \n',Nz);
        fprintf(fid,'type_x = FREE_SLIP \n');
%        fprintf(fid,'type_z = FREE_SLIP \n');
        fprintf(fid,'type_z = NO_SLIP \n');
        fprintf(fid,'min_x = %12.8f\n',0);
        fprintf(fid,'min_z = %12.8f\n',0);
        fprintf(fid,'mapped_grid = false\n');
        
        fprintf(fid,'file_type = MATLAB\n');
        fprintf(fid,'u_file = u.orig\n');
        fprintf(fid,'w_file = w.orig\n');
        fprintf(fid,'rho_file = rho.orig\n');
        
        fprintf(fid,'g = %12.8f \n',g);
        fprintf(fid,'tilt_slope = %12.8f \n',tilting);
	fprintf(fid,'rot_f = %12.8f \n',f);
        fprintf(fid,'rho_0 = %12.8f \n',1000);
        fprintf(fid,'visco = %12.8f \n',nu);
        fprintf(fid,'kappa_rho = %12.8f \n',kappa);
        fprintf(fid,'perturb = %f \n',pert);
        
        fprintf(fid,'final_time = %12.8f\n',final_time);
        fprintf(fid,'plot_interval = %12.8f \n',plot_interval);
        fprintf(fid,'plot_interval_1d = %12.8f \n',plot_interval_1d);
	fprintf(fid,'restart = false \n');
        fclose(fid);
        
        cd ../matlab2spins
    end
end

if ~test, quit, end
