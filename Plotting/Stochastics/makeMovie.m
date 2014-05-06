function makeMovie(stoc,ddft,optsPlot,optsPhys,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        makeMovie1D(stoc,ddft,optsPlot,equilibria)
        
    case 'planar'
        makeMovie1D(stoc,ddft,optsPlot,equilibria)
        
    case 'planar2D'
        makeMoviePlanar2D(stoc,ddft,optsPlot,optsPhys,equilibria)
        
    case 'polar2D'
        %makeMoviePolar2D(stoc,ddft,optsPlot,equilibria)
        makeMoviePlanar2D(stoc,ddft,optsPlot,optsPhys,equilibria)
        
end