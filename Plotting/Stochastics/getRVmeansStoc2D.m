function [meanR,meanV]=getRVmeansStoc2D(x,p,geom,nParticlesS,mS)

    switch geom

        case 'planar2D'

            [meanR,meanV]=getRVmeansStocPlanar2D(x,p,nParticlesS,mS);

        case 'polar2D'

            [meanR,meanV]=getRVmeansStocPolar2D(x,p,nParticlesS,mS);
    end
    
end
