function averagesStruct=getAverages(opts,stocStruct)

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
        rho   = zeros(nBins(1),nBins(2),nSpecies,nPlots);
        v     = zeros(nBins(1),nBins(2),2,nSpecies,nPlots);
        flux  = v;
        boxes = v;
    else
        rho   = zeros(nBins,nSpecies,nPlots);
        v     = rho;
        flux  = v;
        boxes = v;
    end
    
    x = stocStruct.x;
    p = stocStruct.p;

    %----------------------------------------------------------------------
    % Means for each time and species
    %----------------------------------------------------------------------

    if(dim==2)
        [rMean,vMean] = getRVmeansStoc2D(x,p,geom,nParticlesS,mS);
        [~,fluxMean]  = getRFluxMeansStoc2D(x,p,geom,nParticlesS,mS);
    else
        [rMean,vMean] = getRVmeansStoc1D(x,p,geom,dim,nParticlesS,mS);    
        [~,fluxMean]  = getRFluxMeansStoc1D(x,p,geom,dim,nParticlesS,mS);    
    end

    %----------------------------------------------------------------------
    % Distributions for each time and species
    %----------------------------------------------------------------------

    for iPlot = 1: nPlots
        xt=x(:,:,iPlot)'; 
        pt=p(:,:,iPlot)';
        
        [R,P]=getRP(xt,pt,geom,dim);

        if(dim==2)
            if(optsPlot.fixedBins)
                histRange=[optsPlot.rMin(1) optsPlot.rMax(1) optsPlot.rMin(2) optsPlot.rMax(2)];
                [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS,histRange);
            else
                [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS);
            end
            %[nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS);  
            rho(:,:,:,iPlot)   = nR;
            boxes(:,:,:,:,iPlot) = xR;
                        
            for iSpecies=1:nSpecies
                v(:,:,:,iSpecies,iPlot) = meanP(:,:,:,iSpecies)./mS(iSpecies);
                flux(:,:,1,iSpecies,iPlot) = v(:,:,1,iSpecies,iPlot).*rho(:,:,iSpecies,iPlot);
                flux(:,:,2,iSpecies,iPlot) = v(:,:,2,iSpecies,iPlot).*rho(:,:,iSpecies,iPlot);
            end
        else
            [nR,meanP,xR]= RPhist(R,P,nBins,nParticlesS);   
            rho(:,:,iPlot)   = nR;
            boxes(:,:,iPlot) = xR;
            mSrep=repmat(mS',size(meanP,1),1);                
            v(:,:,iPlot)=meanP./mSrep;
            flux(:,:,iPlot) = v(:,:,iPlot).*rho(:,:,iPlot);
        end
    end

    averagesStruct.rMean     = rMean;
    averagesStruct.vMean     = vMean;
    averagesStruct.fluxMean  = fluxMean;
    averagesStruct.rho       = rho;
    averagesStruct.v         = v;
    averagesStruct.flux      = flux;
    averagesStruct.boxes     = boxes;
    averagesStruct.plotTimes = plotTimes;


end