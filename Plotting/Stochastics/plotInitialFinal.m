function plotInitialFinal(stoc,ddft,optsPlot,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        plotInitialFinal1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar'
        plotInitialFinal1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar2D'
        plotInitialFinal2D(stoc,ddft,optsPlot,equilibria);
        
    case 'polar2D'
        plotInitialFinal2D(stoc,ddft,optsPlot,equilibria);
        
end