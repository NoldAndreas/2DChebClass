function plotMeans(stoc,ddft,optsPlot,equilibria)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        plotMeans1D(stoc,ddft,optsPlot,equilibria)
        
    case 'planar'
        plotMeans1D(stoc,ddft,optsPlot,equilibria)
        
    case 'planar2D'
        plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria)
        
    case 'polar2D'
        %plotMeansPolar2D(stoc,ddft,optsPlot,equilibria,pdfFile)
        plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria)
        
end

