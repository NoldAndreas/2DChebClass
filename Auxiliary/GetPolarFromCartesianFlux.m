function fluxPolar = GetPolarFromCartesianFlux(fluxCart,theta)
    ux  = fluxCart(1:end/2,:);
    uy  = fluxCart(1+end/2:end,:);

    [u,v] = GetPolarFromCartesian(ux,uy,theta);
    fluxPolar = [u;v];

end