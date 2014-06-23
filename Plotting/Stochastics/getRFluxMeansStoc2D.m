function [meanR,meanFlux]=getRFluxMeansStoc2D(x,p,geom,nParticlesS,mS)
    
    switch geom

        case 'planar2D'

            [meanR,meanFlux]=getRFluxMeansStocPlanar2D(x,p,nParticlesS,mS);

        case 'polar2D'

            [meanR,meanFlux]=getRFluxMeansStocPolar2D(x,p,nParticlesS,mS);
    end
    
end