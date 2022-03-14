%% driver_icecoveredwaters.m
% Driver script for ICECOVEREDWATERS SPINS case
% icecoveredwaters.txt' contains list of cases to run
% with different parameters. This script makes the corresponding folders
% containing spins.conf files and the executables.
% Other m-files required: none

clearvars

test = true;                       % set to false to write data to disk

% Case independent parameters
% Spatial parameters
Lx = 10.0;
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
type_z = 'NO_SLIP';
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
    amp = 0.01; frq = 2;

% Temporal Parameters
plot_interval = 1;

% Filter Parameters
f_cutoff = 0.6;
f_order = 2.0;
f_strength = 20.0;

% Secondary diagnostics
compute_stresses_bottom = true;
compute_stresses_top = true;
% Case dependent parameters
par = readtable('icecoveredwaters.txt');     % this may not work for R2017b or later

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
    ice_length = par.ice_length(numcase);
    ice_thickness = par.ice_thickness(numcase);
    ice_trans = par.ice_trans(numcase);
    
    x = Lx*(0.5:Nx-0.5)/Nx;     % free-slip in x
    if strcmpi(type_z, 'free_slip')
        z = Lz*(0.5:Nz-0.5)/Nz;     % free-slip in z
    elseif strcmpi(type_z, 'no_slip')
        N = Nz - 1;
        ztmp = cos(pi*(0:N)/N)';       % no-slip in z: size(z1d)==Nz+1
        z = flipud(0.5*Lz*(ztmp+1));    % so that z is in increasing order
    else
        error([type_z, ' Expansion in Z not configured'])
    end
    
    [xx,~] = meshgrid(x,z);	% 2D grid in Matlab's meshgrid format
    ii = repmat(0:Nz-1, Nx,1)';

    % set up Chebyshev in vertical
    zz = -cos(ii.*pi./(Nz-1));
    % Compute zgrid / topography
    %topo = hill_slope/2 *(hill_length + d*(log(2*cosh((xx - a1)/d)) - log(2*cosh((xx-a2)/d))));
    ice_step = ((0.5*tanh((xx - (Lx-ice_length))/ice_trans)) + 0.5);
    roughness = amp*sin(frq*pi*xx).*ice_step;
    
    topo = -(ice_thickness/2 + (ice_thickness*0.5*tanh((xx - (Lx-ice_length))/ice_trans))) + roughness;
    
    zg = min_z + 0.5*Lz*(1+zz) + 0.5*(1+zz).*topo;

    % Compute stratification
    
    rho =  -0.5*delta_rho*tanh((zg-pyc_loc)/h_halfwidth);
    rho = rho.*0.5.*(1.0+tanh((xx-L_adj)/delta_x)); % Clears the region behind the gate
    rho = rho + 0.5*(1.0-tanh((xx-L_adj)/delta_x))...
        .*(-0.5*delta_rho*tanh((zg-pyc_adj_loc)/h_halfwidth)); % Adds in density around the gate    
    
    % Compute initial velocity
    u = zeros(size(xx));
    w = zeros(size(xx));
    
    if test     % plot initial density configuration
        figure(numcase); subplot(3, 1, 1);
        pcolor(xx,zg,rho), shading flat, cmocean('dense'); title('$\rho$', 'interpreter', 'latex');
        subplot(3, 1, 2); pcolor(xx,zg,u), shading flat, cmocean('balance'); title('u')
        subplot(3, 1, 3); pcolor(xx,zg,w), shading flat, cmocean('balance'); title('w')
        
    else        % write data to disk
        mkdir(['../' casename]), cd(['../' casename])
        
        dspins = permute(rho,[2,1]);
        uspins = permute(u,[2,1]);
        wspins = permute(w,[2,1]);
        
        fid = fopen('rho.orig','wb'); fwrite(fid,dspins,'double'); fclose(fid);
        fid = fopen('u.orig','wb'); fwrite(fid,uspins,'double'); fclose(fid);
        fid = fopen('w.orig','wb'); fwrite(fid,wspins,'double'); fclose(fid);
        
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
        
        fprintf(fid,'file_type = MATLAB\n');
        fprintf(fid,'u_file = u.orig\n');
        fprintf(fid,'w_file = w.orig\n');
        fprintf(fid,'rho_file = rho.orig\n');
        
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
        fprintf(fid,'ice_length = %4.3f \n', ice_length);
        fprintf(fid,'ice_thickness = %5.4f \n', ice_thickness);
        fprintf(fid,'ice_trans = %4.3f \n', ice_trans);
        %fprintf(fid,'hill_end_dist = %6.2f \n', hill_end_dist);
        
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
        if compute_stresses_top
            fprintf(fid, 'compute_stresses_top = true');
        end
        fclose(fid);
        
        % Copy executable files across to directory
        if isunix
            copyfile([cdir, '/wave_reader.x'], '.');
            copyfile([cdir, '/derivatives.x'], '.');
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
    fprintf(fid, char(par.casename(i)));
	fprintf(fid, char('\n'));
end
