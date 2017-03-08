function JCP_639_2017_Fig2_GaussianInterpolation    

    AddPaths('JCP_639_2017');    
    close all;

    geom     = struct('N',50,'L',2);
    PlotArea = struct('N',2000,'yMin',-25,'yMax',25);    
    infLine  = InfSpectralLine(geom);            
    infLine.ComputeAll(PlotArea);  

    % **** Plotting **** 
    figure('Position',[0 0 450 200],'color','white');          
    subplot(1,3,1); Plot(geom.N-20,infLine);    
    subplot(1,3,2); Plot(geom.N-6,infLine);    
    subplot(1,3,3); Plot(geom.N-3,infLine);
    SaveFigure('JCP_639_2017_Fig2');
        
    function g = Gaussian(x)
        mu    = 0;
        sigma = 1; 
        g = 1/sqrt(2*pi)/sigma * exp( -(x-mu).^2/(2*sigma^2) );
    end

    function Plot(k,infLine)        
        y  = infLine.Pts.y; 
        infLine.plot(Gaussian(y-y(k)));
        pbaspect([1 1 1]);
        xlim([-25,25]);  ylim([-0.1,0.5]);
        xlabel('$y$','interpreter','latex');
        title(['$g(y-y_{',num2str(k),'})$'],'interpreter','latex');        
    end

end