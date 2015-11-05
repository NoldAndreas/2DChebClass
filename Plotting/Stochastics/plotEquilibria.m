function outputFile = plotEquilibria(stoc,ddft,optsPlotGIF,optsPhys,equilibria)

geom=optsPlotGIF(1).geom;

switch geom         % for different geometries

    case 'spherical'
        outputFile = plotEquilibria1D(stoc,ddft,optsPlotGIF,equilibria);
        %disp('Not implemented for 1D')
        
    case 'planar'
        outputFile = plotEquilibria1D(stoc,ddft,optsPlot,equilibria);
        %disp('Not implemented for 1D')
        
    case 'planar2D'
        outputFile = plotEquilibria2D(stoc,ddft,optsPlotGIF,optsPhys,equilibria);
        
    case 'polar2D'
        outputFile = plotEquilibria2D(stoc,ddft,optsPlotGIF,optsPhys,equilibria);
        
end