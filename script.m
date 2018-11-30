covar.model = 'k-bessel'; % type of covariance function, see functions/kriginginitiaite.m for option
covar.range0 = [10 1000]; % range of covariance [y x]
covar.azimuth = [0]; % orientation of the covariance
covar.c0 = 1; % variance
covar.alpha = 5; % parameter of covariance function (facult)

sim.s=[100 100]; % simulation size 
sim.n=1; % number of simulation

tic
res = simField(sim,covar);
toc

for k=1:sim.n
    imagesc(res{k})
    pause(0.5)
end
