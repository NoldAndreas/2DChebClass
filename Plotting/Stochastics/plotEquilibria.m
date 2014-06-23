function outputFile = plotEquilibria(stoc,ddft,optsPlotGIF,equilibria)

geom=optsPlotGIF(1).geom;

switch geom         % for different geometries

    case 'spherical'
        %plotEquilibria1D(stoc,ddft,optsPlotGIF,equilibria);
        disp('Not implemented for 1D')
        
    case 'planar'
        %plotEquilibria1D(stoc,ddft,optsPlotGIF,equilibria);
        disp('Not implemented for 1D')
        
    case 'planar2D'
        outputFile = plotEquilibria2D(stoc,ddft,optsPlotGIF,equilibria);
        
    case 'polar2D'
        outputFile = plotEquilibria2D(stoc,ddft,optsPlotGIF,equilibria);
        
end