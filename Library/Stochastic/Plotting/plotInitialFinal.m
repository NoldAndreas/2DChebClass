function outputFile = plotInitialFinal(stoc,ddft,optsPlot,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        outputFile = plotInitialFinal1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar'
        outputFile = plotInitialFinal1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar2D'
        outputFile = plotInitialFinal2D(stoc,ddft,optsPlot,equilibria);
        %outputFile = plotInitialFinalError2D(stoc,ddft,optsPlot,equilibria);
        
    case 'polar2D'
        outputFile = plotInitialFinal2D(stoc,ddft,optsPlot,equilibria);
        
end