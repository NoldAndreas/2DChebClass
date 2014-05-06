function plotMeans(stoc,ddft,optsPlot,equilibria,pdfFile)

geom=optsPlot.geom;

switch geom         % for different geometries

    case 'spherical'
        plotMeans1D(stoc,ddft,optsPlot,equilibria,pdfFile)
        
    case 'planar'
        plotMeans1D(stoc,ddft,optsPlot,equilibria,pdfFile)
        
    case 'planar2D'
        plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria,pdfFile)
        
    case 'polar2D'
        %plotMeansPolar2D(stoc,ddft,optsPlot,equilibria,pdfFile)
        plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria,pdfFile)
        
end

