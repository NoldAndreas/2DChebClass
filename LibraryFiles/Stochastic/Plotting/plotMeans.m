function outputFile = plotMeans(stoc,ddft,optsPlot,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        outputFile = plotMeans1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar'
        outputFile = plotMeans1D(stoc,ddft,optsPlot,equilibria);
        
    case 'planar2D'
          
        outputFile = plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria);
        
    case 'polar2D'
        %plotMeansPolar2D(stoc,ddft,optsPlot,equilibria,pdfFile)
        outputFile = plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria);
        
end

