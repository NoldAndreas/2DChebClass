function IntError = IntegrationNL(N,L)

if(nargin==0)
    N = 50;
    L = 2;
end

PhysArea.N = N;
PhysArea.L = L;

aLine = InfSpectralLineSpherical(PhysArea);

PlotArea = struct('N',200,'yMin',0,'yMax',10);

[Pts,~,Int,~,~] = aLine.ComputeAll(PlotArea); 

y = Pts.y;

sigma = 1;
a = 1;

u0Full = (sigma*sqrt(2*pi))^(-3)*exp(-y.^2/(2*sigma^2));

IntError = abs(L1norm(u0Full) - 1);

    function Norm = L1norm(u)
        integrand = y.^2.*abs(u);
        integrand(~isfinite(y)) = 0;
        Norm = 0.5*4*pi*Int*(integrand);
    end

end

