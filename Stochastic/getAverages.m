function averagesStruct=getAverages(opts,stocStruct)

%function equilibria=getEquilibria(xEq,pEq,optsPlot)
%equilibria=getEquilibria(xEq,optsPlot)
%   calculates equilibrium distributions and means from stochastic sampling
%
% INPUTS: 
%  xEq    -- nSamples x (dim*nParticles) matrix of equilibrium positions
%  optsPlot  -- structure of size [1 1] containing
%             plotTimes       (vector of plot times)
%             geom            (geometry; either 'planar' or 'spherical')
%             dim             (dimension)
%             nBins           (number of bins for histogramming)
%             plotDensity     (true/false whether to plot the
%                                 density; false->distruibution)
%
% OUTPUTS:
%  equilbria -- structure of size [1 1] containing
%             REq             (vector of rho at positions xR)
%             vEq             (vector of v at positions xR)
%             pEq             (vector of p at positions xR)
%             xR              (vector positions, centres of histogram bins)
%             rMeanEq         (vector of mean positions, 
%                               same length as plotTimes)
%             vMeanEq         (vector of mean velocities, 
%                               same length as plotTimes)
%             plotTimes       (vector of plot times, same as input)

    %optsPhys = opts.optsPhys;
    %optsStoc = opts.optsStoc;
    optsPlot = opts.optsPlot;
        
    % determine the relevant parameters
    geom=optsPlot.geom;
    nBins=optsPlot.nBins;
    dim=optsPlot.dim;
    plotTimes=optsPlot.plotTimes;    
    nParticlesS=optsPlot.nParticlesS;
    mS=optsPlot.mS;
    
    nSpecies = length(nParticlesS);
    nPlots   = length(plotTimes);
    
    if(dim==2)
        rho   = zeros(nBins(1),nBins(2),2,nSpecies,nPlots);
        v     = rho;
        boxes = rho;
    else
        rho   = zeros(nBins,nSpecies,nPlots);
        v     = rho;
        boxes = rho;
    end
    
    x = stocStruct.x;
    p = stocStruct.p;

    %----------------------------------------------------------------------
    % Means for each time and species
    %----------------------------------------------------------------------

    if(dim==2)
        [rMean,vMean]=getRVmeansStoc2D(x,p,geom,nParticlesS,mS);
    else
        [rMean,vMean]=getRVmeansStoc1D(x,p,geom,dim,nParticlesS,mS);    
    end

    %----------------------------------------------------------------------
    % Distributions for each time and species
    %----------------------------------------------------------------------

    for iPlot = 1: nPlots
        xt=x(:,:,iPlot)'; 
        pt=p(:,:,iPlot)';

        [R,P]=getRP(xt,pt,geom,dim);

        if(dim==2)
            [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS);   
            rho(:,:,:,:,iPlot)   = nR;
            boxes(:,:,:,:,iPlot) = xR;
            for iSpecies=1:nSpecies
                v(:,:,:,iSpecies,iPlot)=meanP(:,:,:,iSpecies)./mS(iSpecies);
            end
        else
            [nR,meanP,xR]= RPhist(R,P,nBins,nParticlesS);   
            rho(:,:,iPlot)   = nR;
            boxes(:,:,iPlot) = xR;
            mSrep=repmat(mS',size(meanP,1),1);                
            v(:,:,iPlot)=meanP./mSrep;
        end
    end

    averagesStruct.rMean     = rMean;
    averagesStruct.vMean     = vMean;
    averagesStruct.rho       = rho;
    averagesStruct.v         = v;
    averagesStruct.boxes     = boxes;
    averagesStruct.plotTimes = plotTimes;


end