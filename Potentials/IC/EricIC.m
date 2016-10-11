function rho = EricIC(r,optsPhys)

    N = optsPhys.nParticles;

    rho = N*(2*pi)^(-1/2)*exp(-r.^2/2);

end