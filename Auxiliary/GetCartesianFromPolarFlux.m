function fluxCart = GetCartesianFromPolarFlux(fluxPolar,theta)
    ur  = fluxPolar(1:end/2,:);
    uth = fluxPolar(1+end/2:end,:);

    [u,v] = GetCartesianFromPolar(ur,uth,theta);
    fluxCart = [u;v];

end