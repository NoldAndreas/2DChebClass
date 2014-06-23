function equilibria=getEquilibria(opts,xpStruct)

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

    xEq = xpStruct.xEq;
    pEq = xpStruct.pEq;
    
    
    % determine the relevant parameters
    geom=optsPlot.geom;
    nBins=optsPlot.nBins;
    dim=optsPlot.dim;
    plotTimes=optsPlot.plotTimes;
    
    nParticlesS=optsPlot.nParticlesS;
    
    % species masses
    mS=optsPlot.mS;
    
    %----------------------------------------------------------------------
    % Means for each time and species
    %----------------------------------------------------------------------
    
    % get the mean 1D positions and momenta for each species
        
    if(optsPlot.dim==2)
        [rMeanEq,vMeanEq]=getRVmeansStoc2D(xEq,pEq,geom,nParticlesS,mS);
    else
        [rMeanEq,vMeanEq]=getRVmeansStoc1D(xEq,pEq,geom,dim,nParticlesS,mS);    
    end
    
    % form matrices to be plotted in time
    % rMeanEq is of size 1 x nSpecies x Dim
    % rMeanEq should be of size plotTimes x nSpecies x Dim
    
    rMeanEq=repmat(rMeanEq,[length(plotTimes) 1 1]);
    vMeanEq=repmat(vMeanEq,[length(plotTimes) 1 1]);
 
    %----------------------------------------------------------------------
    % Distributions for each species
    %----------------------------------------------------------------------
    
    % getRP takes matrices of size (dim*nParticles) x nSamples
    xEq=xEq';
    pEq=pEq';

    % get the R and P data
    [R,P]=getRP(xEq,pEq,geom,dim);
    
    if(optsPlot.dim==2)
        histRange=[optsPlot.rMin(1) optsPlot.rMax(1) optsPlot.rMin(2) optsPlot.rMax(2)];
        [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS,histRange);
        meanV=meanP;
        for iSpecies=1:length(mS)
            meanV(:,:,:,iSpecies)=meanP(:,:,:,iSpecies)./mS(iSpecies);
        end
    else
        [nR,meanP,xR]= RPhist(R,P,nBins,nParticlesS);   
        mS=repmat(mS',size(meanP,1),1);
        meanV=meanP./mS;
    end    


    % form output structure
    equilibria.REq=nR;
    equilibria.vEq=meanV;
    equilibria.pEq=meanP;
    equilibria.xEq=xR;
    equilibria.rMeanEq=rMeanEq;
    equilibria.vMeanEq=vMeanEq;
    equilibria.plotTimes=plotTimes;

end