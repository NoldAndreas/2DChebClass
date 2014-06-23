function outputFile = makeMovie(stoc,ddft,optsPlot,optsPhys,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        outputFile = makeMovie1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar'
        outputFile = makeMovie1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar2D'
        outputFile = makeMoviePlanar2D(stoc,ddft,optsPlot,optsPhys,equilibria);
        
    case 'polar2D'
        %makeMoviePolar2D(stoc,ddft,optsPlot,equilibria)
        outputFile = makeMoviePlanar2D(stoc,ddft,optsPlot,optsPhys,equilibria);
        
end