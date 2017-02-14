function JCP_639_2017_Fig1_IntegrationNLPlots

    AddPaths('JCP_639_2017');
    close all

    Nvals  = 10:10:200;
    Lvals  = 1:14;
    IntErr = zeros(length(Nvals),length(Lvals));

    for iN = 1:length(Nvals)
        for iL = 1:length(Lvals)
            IntErr(iN,iL) = IntegrationNL(Nvals(iN),Lvals(iL));
        end
    end

    %**** Plotting ****        
    hIntErr    = figure('Position',[0 0 1000 800]);
    imagesc(log10(IntErr));
    colormap(hIntErr,gray)

    set(gca,'YTick',1:length(Nvals));
    set(gca,'YTickLabel',Nvals);
    set(gca,'XTick',1:length(Lvals));
    set(gca,'XTickLabel',Lvals);
    xlabel('$L$','interpreter','latex');
    ylabel('$N$','interpreter','latex');
    SaveFigure('JCP_639_2017_Fig1_IntegrationNLPlots');

    function IntError = IntegrationNL(N,L)

        PhysArea   = struct('N',N,'L',L);
        PlotArea   = struct('N',200,'yMin',0,'yMax',10);

        aLine      = InfSpectralLineSpherical(PhysArea);       
        [Pts,~,Int,~,~] = aLine.ComputeAll(PlotArea); 

        y     = Pts.y;
        sigma = 1;

        u0Full = (sigma*sqrt(2*pi))^(-3)*exp(-y.^2/(2*sigma^2));
        IntError = abs(L1norm(u0Full) - 1);

        function Norm = L1norm(u)
            integrand = y.^2.*abs(u);
            integrand(~isfinite(y)) = 0;
            Norm = 0.5*4*pi*Int*(integrand);
        end

    end

end