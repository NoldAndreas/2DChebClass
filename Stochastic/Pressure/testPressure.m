clear all;

optsPhys.nParticlesS = 10;
optsPhys.nParticles = optsPhys.nParticlesS;
optsPhys.dim = 1;
optsPhys.kBT = 1;
optsPhys.t = 0;
optsPhys.m = 1;
optsPhys.geom       = 'full';
optsPhys.type       = 'stoc';

optsPhys.sigma = 0;
optsPhys.L = 5;


optsStoc.nSamples = 20000;
optsStoc.thin = 1;
optsStoc.burnin = 10000;

optsPhys.V1DV1      = 'MB';
optsPhys.V2DV2      = 'free';
optsStoc.initialGuess = 'MBIG';

opts.optsPhys = optsPhys;
opts.optsStoc = optsStoc;

P = samplepdf(opts,[]);

%optsPhys.V1DV1      = 'zeroPotential1D';
%optsPhys.V1DV1      = 'V1_SlitStoc';
optsPhys.V1DV1      = 'quadratic1D';
optsPhys.Vm = 1;
%optsPhys.V2DV2      = 'Gaussian';
optsPhys.V2DV2      = 'free';
optsPhys.epsilon = 1;
optsPhys.alpha   = 1;

optsPhys.rMin = -10;
optsPhys.rMax = 10;


optsStoc.initialGuess = 'MBIG';

opts.optsPhys = optsPhys;

R = samplepdf(opts,[]);

P = P';
R = R';

nBins = 20;

[PK,PV,binMids,nR,meanP] = getPressure1D(R,P,nBins,optsPhys)
