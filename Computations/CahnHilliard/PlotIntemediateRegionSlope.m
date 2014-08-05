function PlotIntemediateRegionSlope

    x  = (0.2:0.01:1)';
    Ca = 0.03;
    thetam = 0.1;

    z =  Ca*log(x) + 0.3;

    plot(x,IntermediateStokesSlope(z,1,thetam));

end