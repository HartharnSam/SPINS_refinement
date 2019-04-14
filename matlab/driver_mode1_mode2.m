%% driver_mode1_mode.m
% Driver script for mode1/mode2 tutorial from the fluids wiki.
% Parameter file 'mode1_mode2.txt' contains list of cases to run
% with different parameters. This script makes the corresponding
% spins.conf files and the initial u,v,w,rho fields for the model.

clearvars

test = false;                       % set to false to write data to disk

% Case independent parameters
min_x=0.0;
min_y=0.0;
min_z=0.0;

type_x='Free_SLIP';
type_y='FOURIER';
type_z='FREE_SLIP';

g=9.81;
rot_f=0.0;
rho_0=1000.0;
visco=1.e-6;
kappa_rho=2.e-7;
kappa_tracer=2.e-7;

delta_rho=0.02;
h_halfwidth=0.01;

pert = 0.0;
%% Get rid of this block and just put it into SPINS.conf?
h_halfwidth=0.01;
delta_x=0.04;
dye1_loc=2.0;
dye2_loc=5.0;
dye1_width=0.2;
dye2_width=0.2;
dye_halfwidth=0.005;
enable_tracer=true;

final_time=60;
plot_interval=1;

% Case dependent parameters
par = readtable('mode1_mode2.txt');     % this may not work for R2017b or later

cdir=pwd; % Remember directory script is run from
%% Loop through different cases
for numcase = 1:height(par)
    casename = char(par.casename(numcase));

    % Grid
    Lx = par.Lx(numcase);
    Lz = par.Lz(numcase);
    Nx = par.Nx(numcase);
    Nz = par.Nz(numcase);
    
    pyc_loc=par.pyc_loc(numcase);
    H1=par.H1(numcase);
    H2=par.H2(numcase);
    L1=par.L1(numcase);
    L2=par.L2(numcase);

    x = Lx*(0.5:Nx-0.5)/Nx;     % free-slip in x
    z = Lz*(0.5:Nz-0.5)/Nz;     % free-slip in z
%    [~,ztmp] = cheb(Nz);            % no-slip in z: size(z1d)==Nz+1
%    z = flipud(0.5*Lz*(ztmp+1));    % so that z is in increasing order
    [xx,zz] = meshgrid(x,z);	% 2D grid in Matlab's meshgrid format
    
    % Stratification
    rho1 = -0.5*delta_rho*tanh( (zz - (pyc_loc - H1))/h_halfwidth);
    rho2 = -0.5*delta_rho*tanh((zz-pyc_loc)/h_halfwidth);
    rho3 = -0.25*delta_rho* ( tanh((zz-pyc_loc + 0.5*H2)/h_halfwidth) + tanh((zz-pyc_loc - 0.5*H2)/h_halfwidth));

    filt1=0.5*(1 - tanh((xx-L1)/delta_x));
    filt2=0.5*(tanh((xx-L1)/delta_x) - tanh((xx-(Lx-L2))/delta_x));
    filt3=0.5*(1+tanh((xx-(Lx-L2))/delta_x));

    rho = filt1.*rho1 + filt2.*rho2 + filt3.*rho3;

    
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
        fprintf(fid,'Ly = 0.30 \n');
        fprintf(fid,'Lz = %12.8f \n',Lz);
        fprintf(fid,'Nx = %d \n',Nx);
        fprintf(fid,'Ny = 1 \n');
        fprintf(fid,'Nz = %d \n',Nz);

        fprintf(fid,'min_x = 0 \n');
        fprintf(fid,'min_y = 0 \n');
        fprintf(fid,'min_z = 0 \n');

        fprintf(fid,'type_x = FREE_SLIP \n');
        fprintf(fid,'type_y = FOURIER \n');
        fprintf(fid,'type_z = FREE_SLIP \n');
%        fprintf(fid,'type_z = NO_SLIP \n');
        fprintf(fid,'mapped_grid = false\n');
        
        fprintf(fid,'file_type = MATLAB\n');
        fprintf(fid,'u_file = u.orig\n');
        fprintf(fid,'w_file = w.orig\n');
        fprintf(fid,'rho_file = rho.orig\n');
        
        fprintf(fid,'g = %12.8f \n',g);
	    fprintf(fid,'rot_f = %12.8f \n',rot_f);
        fprintf(fid,'rho_0 = %12.8f \n',rho_0);
        fprintf(fid,'visco = %12.8f \n',visco);
        fprintf(fid,'kappa_rho = %12.8f \n',kappa_rho);
        fprintf(fid,'kappa_tracer = 2e-7 \n');
        fprintf(fid,'perturb = %f \n',pert);

        fprintf(fid,'delta_rho = %12.8f \n',delta_rho);
        fprintf(fid,'pyc_loc = %12.8f \n',pyc_loc);
        fprintf(fid,'H1=%12.8f \n',H1);
        fprintf(fid,'H2=%12.8f \n',H2);
        fprintf(fid,'L1=%12.8f \n',L1);
        fprintf(fid,'L2=%12.8f \n',L2);
        fprintf(fid,'h_halfwidth=%12.8f \n',h_halfwidth);
        fprintf(fid,'delta_x=%12.8f \n',delta_x);
%        fprintf(fid,'dye1_loc=%12.8f \n',dye1_loc);
%        fprintf(fid,'dye2_loc=%12.8f \n',dye2_loc);
%        fprintf(fid,'dye1_width=%12.8f \n',dye1_width);
%        fprintf(fid,'dye2_width=%12.8f \n',dye2_width);
%        fprintf(fid,'dye_halfwidth=0.005 \n');
        fprintf(fid,'enable_tracer = false \n');

%        fprintf(fid,'hill_height = 0.0 \n');
%        fprintf(fid,'hill_centre = 4.8 \n');
%        fprintf(fid,'hill_width = 0.10 \n');
        
        fprintf(fid,'final_time = %12.8f\n',final_time);
        fprintf(fid,'plot_interval = %12.8f \n',plot_interval);
        fprintf(fid,'plot_interval_1d = %12.8f \n',plot_interval);
	    
        fprintf(fid,'restart = false \n');
        fprintf(fid,'restart_time = 0.0 \n');
        fprintf(fid,'restart_sequence=0 \n');
        fprintf(fid,'restart_from_dump = false \n');
        fprintf(fid,'compute_time = -1 \n');

        fprintf(fid,'f_cutoff = 0.6 \n');
        fprintf(fid,'f_order = 2.0 \n');
        fprintf(fid,'f_strength = 20.0 \n');
        fclose(fid);
        
        cd(cdir);
    end
end



if ~test, quit, end
