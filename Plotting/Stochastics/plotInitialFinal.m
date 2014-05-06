function plotInitialFinal(stoc,ddft,optsPlotGIF,equilibria)

geom=optsPlotGIF.geom;

switch geom         % for different geometries

    case 'spherical'
        plotInitialFinal1D(stoc,ddft,optsPlotGIF,equilibria);
        
    case 'planar'
        plotInitialFinal1D(stoc,ddft,optsPlotGIF,equilibria);
        
    case 'planar2D'
        plotInitialFinal2D(stoc,ddft,optsPlotGIF,equilibria);
        
    case 'polar2D'
        plotInitialFinal2D(stoc,ddft,optsPlotGIF,equilibria);
        
end