function plotEquilibria(stoc,ddft,optsPlotGIF,equilibria,pdfFile)

geom=optsPlotGIF(1).geom;

switch geom         % for different geometries

    case 'spherical'
        %plotEquilibria1D(stoc,ddft,optsPlotGIF,equilibria,pdfFile);
        disp('Not implemented for 1D')
        
    case 'planar'
        %plotEquilibria1D(stoc,ddft,optsPlotGIF,equilibria,pdfFile);
        disp('Not implemented for 1D')
        
    case 'planar2D'
        plotEquilibria2D(stoc,ddft,optsPlotGIF,equilibria,pdfFile);
        
    case 'polar2D'
        plotEquilibria2D(stoc,ddft,optsPlotGIF,equilibria,pdfFile);
        
end