function Simulation1DInterpolation
    
    N = 20;
    L = 1;

    x       = ClenCurtFlip(N);
    y       = QuotientMap(x,L,1,Inf);
    xInterp = -1:0.01:0.9;
    yInterp = QuotientMap(xInterp,L,1,Inf);
    IP      = barychebevalMatrix(x,xInterp');    
    z       = BarkerHenderson_2D(y)';
    
    plot(y,z,'o','MarkerFaceColor','g'); hold on;
    plot(yInterp,IP*z);
    xlim([0 10]);
        
end