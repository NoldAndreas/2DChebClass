function [maxRelErr,t,DError,D2Error,DDError,IntError] = SphericalHeatEquation(N,L)

if(nargin==0)
    N = 50;
    L = 2;
end

tMax = 2;

plotTimes = 0:tMax/20:tMax;

PhysArea.N = N;
PhysArea.L = L;

aLine = InfSpectralLineSpherical(PhysArea);

PlotArea = struct('N',200,'yMin',0,'yMax',10);

tic
[Pts,Diff,Int,Ind,Interp] = aLine.ComputeAll(PlotArea); 

Dy = Diff.Dy;
DDy = Diff.DDy;
y = Pts.y;
yCut = cut(y);

sigma = 1;
a = 1;

u0Full = (sigma*sqrt(2*pi))^(-3)*exp(-y.^2/(2*sigma^2));
u0 = cut(u0Full);

DError = L1norm(Dy*u0Full-DuExact(y,0))/L1norm(DuExact(y,0));

D2Error = L1norm(DDy*u0Full-DDuExact(y,0))/L1norm(DDuExact(y,0));

DDError = L1norm(Dy*(Dy*u0Full)-DDuExact(y,0))/L1norm(DDuExact(y,0));

IntError = abs(L1norm(u0Full) - 1);

opts = odeset('RelTol',5*10^-14,'AbsTol',5*10^-14);

[~,U_t] =  ode15s(@du_dt,plotTimes,u0,opts);
t=toc


relErr = zeros(size(plotTimes));

subplot(1,2,1);
hold off
cla
subplot(1,2,2);
hold off
cla

for iPlot = 1:length(plotTimes);
    subplot(1,2,1)
    ut = mirror(U_t(iPlot,:).');
    plot(Interp.pts,Interp.InterPol*ut,'b');
    hold on
    plot(yCut,U_t(iPlot,:),'go');
    plot(Interp.pts,uExact(Interp.pts,plotTimes(iPlot)),'-r');
    ylim([0,0.1]);
    xlim([0,max(Interp.pts)]);
    
    uEx = uExact(y,plotTimes(iPlot));
    
    subplot(1,2,2);
    %m = L1norm(ut);

    if(iPlot == 1)
        relErr(iPlot) = 0;
    else
        relErr(iPlot) = L1norm(ut-uEx);
    end
    hold on
    plot(plotTimes(iPlot),relErr(iPlot),'x');
    
end


maxRelErr = max(relErr);

    function dudt = du_dt(t,u)
        u = mirror(u);
        dudt = 2*(Dy*u)./y + DDy*u;
        dudt = a^2*dudt;
        %dudt(Ind.infinite,:) = u(Ind.infinite,:) - u0Full(Ind.infinite,:);
        dudt(Ind.infinite,:) = 0;
        dudt = cut(dudt);
        dudt = dudt(:);
    end

    function xOut = mirror(x)
        % flip for negative spatial part
        xOut = [flipdim(x,1); x];
    end

    function xOut = cut(x)
        % remove negative spatial part
        xOut = x(end/2+1:end,:);
    end

    function u = uExact(r,t)
        u = (2*pi*(2*t*a^2 + sigma^2))^(-3/2) ...
            * exp( -r.^2 / (2*(2*t*a^2 + sigma^2)));
    end

    function Du = DuExact(r,t)
        Du = -r/(2*t*a^2 + sigma^2).*uExact(r,t);
    end

    function DDu = DDuExact(r,t)
        DDu = (r.^2-2*t*a^2-sigma^2)/(2*t*a^2 + sigma^2)^2.*uExact(r,t);
    end


    function Norm = L1norm(u)
        integrand = y.^2.*abs(u);
        integrand(~isfinite(y)) = 0;
        Norm = 0.5*4*pi*Int*(integrand);
    end

end

