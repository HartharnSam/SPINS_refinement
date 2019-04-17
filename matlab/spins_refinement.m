% SPINS grid refinement for 2D domain with no-slip conditions and a flat
% bottom boundary.

clearvars
nx = 4;           % nx times original resolution in x
nz = 2;           % nz times original resolution in z
ii = 100;         % number of output to interpolate
method = 'spline';      % Matlab built-in interpolation method
                        % 'nearest', 'linear', 'spline' or 'cubic'

%% Read in old spins.conf, grids and fields
params = spins_params();
x = xgrid_reader();
z = zgrid_reader();
den = spins_reader('rho',ii);
u = spins_reader('u',ii);
w = spins_reader('w',ii);

%transpose so the ordering is correct
x = x';
z = z';
u = u';
w = w';
den = den';

Lx = params.Lx;
Lz = params.Lz;
Nx = params.Nx;
Nz = params.Nz;
type_x = params.type_x;
type_z = params.type_z;
dx = x(2,2) - x(1,1);

%% Define new grids
Nx_new = nx*Nx;      % n times resolution in x
Nz_new = nz*Nz;      % n times resolution in z

% Update grids in spins.conf file
params.Nx = Nx_new;
params.Nz = Nz_new;

% Compute new grids
x1d = Lx*(0.5:Nx_new-0.5)/Nx_new;           % periodic/free-slip in x
switch type_z
    case 'FREE_SLIP'
        z1d = Lz*(0.5:Nz_new-0.5)/Nz_new;	% free-slip in z
    case 'NO_SLIP'
        [~,ztmp] = cheb(Nz_new);            % no-slip in z: size(z1d)==Nz+1
        z1d = flipud(0.5*Lz*(ztmp+1));      % so that z1d is in increasing order
end

[x2d,z2d] = meshgrid(x1d,z1d);      % 2D grid in Matlab's meshgrid format

%% Interpolation

% Expand the grid first
x1 = [x(:,1)-2*dx	x(:,1)-dx   x	x(:,end)+dx   x(:,end)+2*dx];
z1 = [z(:,1)        z(:,1)      z   z(:,end)      z(:,end)];
x1 = [x1(1,:);      x1(1,:);    x1; x1(end,:);    x1(end,:)];
z1 = [2*z1(1,:)-z1(3,:);        2*z1(1,:)-z1(2,:);      z1; ...
    2*z1(end,:)-z1(end-1,:);	2*z1(end,:)-z1(end-2,:)];

% Extrapolate beyond the boundaries
switch type_x
    case 'FREE_SLIP'
        % Dirichlet for u, Neumann for v, w and den in x
        u1   = [ -u(:,2)   -u(:,1)   u    -u(:,end)   -u(:,end-1)];
        w1   = [  w(:,2)    w(:,1)   w     w(:,end)    w(:,end-1)];
        den1 = [den(:,2)  den(:,1)  den  den(:,end)  den(:,end-1)];
    case 'FOURIER'
        % Periodic in x-direction
        u1   = [  u(:,(end-1):end)   u     u(:,1:2)];
        w1   = [  w(:,(end-1):end)   w     w(:,1:2)];
        den1 = [den(:,(end-1):end)  den  den(:,1:2)];
end

% Free/no-slip (Dirichlet for w, Neumann for u, v and den) in z
u1   = [  u1(2,:);   u1(1,:);  u1;    u1(end,:);   u1(end-1,:)];
w1   = [ -w1(2,:);  -w1(1,:);  w1;   -w1(end,:);  -w1(end-1,:)];
den1 = [den1(2,:); den1(1,:); den1; den1(end,:); den1(end-1,:)];
%%
% Interpolation
u2   = interp2(x1,z1,u1,x2d,z2d,method);
w2   = interp2(x1,z1,w1,x2d,z2d,method);
rho2 = interp2(x1,z1,den1,x2d,z2d,method);

%% Write data to disk

dspins = permute(rho2,[2,1]);
uspins = permute(u2,[2,1]);
wspins = permute(w2,[2,1]);

% Make a subdirectory and put new stuff there
mkdir('high_res'), cd('high_res')

fid = fopen('rho.orig','wb'); fwrite(fid,dspins,'double'); fclose(fid);
fid = fopen('u.orig','wb'); fwrite(fid,uspins,'double'); fclose(fid);
fid = fopen('w.orig','wb'); fwrite(fid,wspins,'double'); fclose(fid);

%% Write parameters to spins.conf

fid = fopen('spins.conf','wt');

for i=1:numel(fieldnames(params))
    fields=fieldnames(params);
    fieldname=char(fields(i));
    value=string(params.(fieldname));
    fprintf(fid,'%s = %s \n',fieldname,value);
end

fclose(fid);

%quit
