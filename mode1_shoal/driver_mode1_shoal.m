%% driver_mode1_shoal.m
% Driver script for mode1 shoal case.
% Parameter file 'mode1_shoal.txt' contains list of cases to run
% with different parameters. This script makes the corresponding
% spins.conf files and the initial u,v,w,rho fields for the model.

clearvars

test = false;                       % set to false to write data to disk

% Case independent parameters
% Spatial parameters
Lx = 7.0;
Ly = 0.1;
Lz = 0.3;
Nx = 4096; 
Ny = 1;
Nz = 256;
min_x = 0.0;
min_y = 0.0;
min_z = -0.3;
% Expansion types
type_x = 'FREE_SLIP';
type_y = 'FOURIER';
type_z = 'FREE_SLIP';
mapped_grid = 'true';
% Physical parameters
g = 9.81;
rot_f = 0.0;
rho_0 = 1026.0;
visco = 1.e-6;
kappa_rho = 1.e-7;
kappa_tracer = 1.e-7;

pert = 1.e-3;

%% These are specific to this experiment
delta_rho = 0.019;
delta_x = 0.04;
L_adj = 0.3;
dye_thickness = 0.05;
dye_halfwidth = 0.01;

% Temporal Parameters
plot_interval = 1;

% Filter Parameters
f_cutoff = 0.6;
f_order = 2.0;
f_strength = 20.0;

% Topography parameters
hill_end_dist = -0.0;

% Secondary diagnostics
compute_stresses_bottom = true;

% Case dependent parameters
par = readtable('mode1_shoal.txt');     % this may not work for R2017b or later

cdir=pwd; % Remember directory script is run from
%% Loop through different cases
for numcase = 1:height(par)
    casename = char(par.casename(numcase));

    %% Read in case dependent parameters here
    
    % Problem Parameters
    pyc_adj_loc = par.pyc_adj_loc(numcase);
    pyc_loc = par.pyc_loc(numcase);
    h_halfwidth = par.h_halfwidth(numcase);
    final_time = par.final_time(numcase);
    
    % Problem Topography Parameters
    hill_height = par.hill_height(numcase);
    hill_slope = par.hill_slope(numcase);
    hill_trans = par.hill_trans(numcase);

    x = Lx*(0.5:Nx-0.5)/Nx;     % free-slip in x
    if strcmpi(type_z, 'free_slip')
        z = Lz*(0.5:Nz-0.5)/Nz;     % free-slip in z
    elseif strcmpi(type_z, 'no_slip')
        [~,ztmp] = cheb(Nz);            % no-slip in z: size(z1d)==Nz+1
        z = flipud(0.5*Lz*(ztmp+1));    % so that z is in increasing order
    else 
        error([type_z, ' Expansion in Z not configured'])
    end
    [xx,zz] = meshgrid(x,z);	% 2D grid in Matlab's meshgrid format
    
    % Compute stratification
    %TODO

    
    % Compute initial velocity
    u = zeros(size(xx));
    w = zeros(size(xx));

    if test     % plot initial density configuration
        figure(numcase); subplot(3, 1, 1);
        pcolor(xx,zz,rho), shading flat, colormap darkjet
        subplot(3, 1, 2); pcolor(xx,zz,u), shading flat, colormap darkjet
        subplot(3, 1, 3); pcolor(xx,zz,2), shading flat, colormap darkjet

    else        % write data to disk
        mkdir(['../' casename]), cd(['../' casename])
        
        %% Write params to spins.conf. Change what is written according to your experiment.
        fid = fopen('spins.conf','wt');
        fprintf(fid, '## MODE-1 ISW Shoaling Configuration File \n');
        fprintf(fid, 'name = %s', casename);
        fprintf(fid, '\n # Spatial Parameters \n');
        fprintf(fid,'Lx = %6.2f \n', Lx);
        fprintf(fid,'Ly = %6.2f \n', Ly);
        fprintf(fid,'Lz = %6.2f \n', Lz);
        fprintf(fid,'Nx = %d \n', Nx);
        fprintf(fid,'Ny = %d \n', Ny);
        fprintf(fid,'Nz = %d \n', Nz);
        
        fprintf(fid,'min_x = %6.2f \n', min_x);
        fprintf(fid,'min_y = %6.2f \n', min_y);
        fprintf(fid,'min_z = %6.2f \n', min_z);

        fprintf(fid, '\n # Expansion types \n');
        fprintf(fid,'type_x = %s \n', type_x);
        fprintf(fid,'type_y = %s \n', type_y);
        fprintf(fid,'type_z = %s \n', type_z);
        fprintf(fid,'mapped_grid = %s \n', mapped_grid);
                
        fprintf(fid, '\n # Physical Parameters \n');
        fprintf(fid,'g = %12.3f \n', g);
	    fprintf(fid,'rot_f = %12.3f \n', rot_f);
        fprintf(fid,'rho_0 = %12.1f \n', rho_0);
        fprintf(fid,'visco = %.2e \n', visco);
        fprintf(fid,'kappa_rho = %.2e \n', kappa_rho);
        fprintf(fid,'kappa_tracer = %.2e \n', kappa_tracer);
        fprintf(fid,'perturb = %.2e \n', pert);

        fprintf(fid, '\n # Problem Parameters \n');
        fprintf(fid,'delta_rho = %4.3f \n', delta_rho);
        fprintf(fid,'pyc_loc = %4.3f \n', pyc_loc);
        fprintf(fid,'h_halfwidth =%4.3f \n',h_halfwidth);
        fprintf(fid,'pyc_adj_loc =%4.3f \n',pyc_adj_loc);
        fprintf(fid,'h_adj_halfwidth = %4.3f \n', h_halfwidth);
        fprintf(fid,'delta_x = %4.3f \n',delta_x);
        fprintf(fid,'L_adj = %4.3f \n',L_adj);
        fprintf(fid,'dye_thickness=%4.3f \n',dye_thickness);
        fprintf(fid,'dye_halfwidth = %4.3f \n', dye_halfwidth);
        
        fprintf(fid, '\n # Topography Parameters \n');
        fprintf(fid,'hill_height = %4.3f \n', hill_height);
        fprintf(fid,'hill_slope = %5.4f \n', hill_slope);
        fprintf(fid,'hill_trans = %4.3f \n', hill_trans);
        fprintf(fid,'hill_end_dist = %6.2f \n', hill_end_dist);
        
        fprintf(fid, '\n # Temporal Parameters \n');
        fprintf(fid,'final_time = %12.8f\n',final_time);
        fprintf(fid,'plot_interval = %12.8f \n',plot_interval);
        
        fprintf(fid, '\n # Restart Parameters \n');
        fprintf(fid,'restart = false \n');
        fprintf(fid,'restart_time = 0.0 \n');
        fprintf(fid,'restart_sequence=0 \n');
        fprintf(fid,'restart_from_dump = false \n');
        fprintf(fid,'compute_time = -1 \n');

        fprintf(fid, '\n # Filter Parameters \n');
        fprintf(fid,'f_cutoff = %12.8f \n', f_cutoff);
        fprintf(fid,'f_order = %12.8f \n', f_order);
        fprintf(fid,'f_strength = %12.8f \n', f_strength);
        
        fprintf(fid, '\n # Diagnostics \n');
        if compute_stresses_bottom
            fprintf(fid, 'compute_stresses_bottom = true');
        end
        
        fclose(fid);
        
        % Copy executable files across to directory
        if isunix
            copyfile([cdir, './mode1_shoal.x'], '.');
            copyfile([cdir, './derivatives.x'], '.');
        end

        % And make the derivatives config file
        fid = fopen('spins.conf_deriv','wt');
        fprintf(fid, '## derivative Configuration File \n');
        fprintf(fid, '\n # Spatial Parameters \n');
        fprintf(fid,'Lx = %6.2f \n', Lx);
        fprintf(fid,'Ly = %6.2f \n', Ly);
        fprintf(fid,'Lz = %6.2f \n', Lz);
        fprintf(fid,'Nx = %d \n', Nx);
        fprintf(fid,'Ny = %d \n', Ny);
        fprintf(fid,'Nz = %d \n', Nz);
        
        fprintf(fid, '\n # Expansion types \n');
        fprintf(fid,'type_x = %s \n', type_x);
        fprintf(fid,'type_y = %s \n', type_y);
        fprintf(fid,'type_z = %s \n', type_z);
        fprintf(fid,'mapped_grid = %s \n', mapped_grid);
        if Ny == 1
            fprintf(fid,'v_exist = false \n');
        else
            warning('derivatives not properly configured');
        end
        
        fprintf(fid, '\n # Physical Parameters \n');
        fprintf(fid,'rho_0 = %12.1f \n', rho_0);
        fprintf(fid,'visco = %.2e \n', visco);
        
        fprintf(fid, '\n # Derivative Options \n');
        fprintf(fid,'deriv_files = u \n');
        fprintf(fid,'start_sequence = %d \n', 0);
        fprintf(fid,'final_sequence = %d \n', final_time);
        fprintf(fid,'step_sequence = %d \n', 1);
        fprintf(fid,'deriv_x = false \n');
        fprintf(fid,'deriv_y = false \n');
        fprintf(fid,'deriv_z = false \n');
        fprintf(fid,'do_vor_x = false \n');
        fprintf(fid,'do_vor_y = true \n');
        fprintf(fid,'do_vor_z = false \n');
        fprintf(fid,'do_enstrophy = false \n');
        fprintf(fid,'do_dissipation = true \n');
        fprintf(fid,'do_vort_stretch = false \n');
        fprintf(fid,'do_enst_stretch = false \n');

        fclose(fid);
        
        cd(cdir);       

    end
end

fid = fopen('../case_list','wt');
for i = 1:height(par)
   fprintf(fid, char(par.casename(numcase)));
end