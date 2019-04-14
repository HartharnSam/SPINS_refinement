addpath(genpath('~/projects/ctb-mmstastn/tghill/SPINSmatlab'));

home=pwd;

dirnames=["base","mode1","mode2","lamp","gravity"];
for j=1:length(dirnames)
cd(char(dirnames(j)))
mkdir figures;
gdpar=spins_gridparams();
for i=0:10:gdpar.params.noutputs
    spins_plot2d('rho',i,'savefig',true,'clim',[-0.01,0.01]);
end
cd(home)
end
