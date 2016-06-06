function plotFiles = doPlots(stocStruct,DDFTStruct,optsStocFull,optsDDFTFull,optsPlot,optsPlotParticles,optsPhys)

%function doPlots(stocStruct,DDFTStruct,optsPlotGIF,optsPlotParticles,fileStruct,optsStruct,optsPhysGIF)
%doPlots(stocStruct,DDFTStruct,optsPlotGIF,optsPlotParticles,fileStruct,optsStruct)
%  Determines which plots we're doing and does them
%
% INPUTS:
%          stocStruct            -- structure containing stochastic data [see doStoc]
%          DDFTStruct            -- structure containing DDFT data [see doDDFT] 
%          optsPlotGIF           -- structure containing standard plotting
%                                    options [see getOptsPlot]
%          optsPlotParticles     -- structure containing particle plotting
%                                    options [see getOptsPlot]]
%          optsStruct            -- structure containing all other
%                                    information [see preProcess]
%
% OUTPUTS:
%          <none>


% use this if matlab crashes on 3D plotting:
% opengl software

%--------------------------------------------------------------------------
% Set up output file list
%--------------------------------------------------------------------------

plotFiles = {};

anyPlots  = ~isempty(optsPlot);
anyPlotsP = ~isempty(optsPlotParticles);

%--------------------------------------------------------------------------
% Set text size and line width
%--------------------------------------------------------------------------

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',2)

if(~isempty(stocStruct))

    %--------------------------------------------------------------------------
    % get initial equilibrium quantities to add to plots
    %--------------------------------------------------------------------------
    
    optsStoc = optsStocFull(1);

    fixedInitial = optsStoc.fixedInitial;

    optsStoc = rmfield(optsStoc,'loadSamples');
    optsStoc = rmfield(optsStoc,'fixedInitial');
    optsStoc = rmfield(optsStoc,'sampleFinal');
    optsStoc = rmfield(optsStoc,'poolsize');
    optsStoc = rmfield(optsStoc,'doStoc');
    optsStoc = rmfield(optsStoc,'loadStoc');
    optsStoc = rmfield(optsStoc,'saveStoc');
    optsStoc = rmfield(optsStoc,'HI');
    optsStoc = rmfield(optsStoc,'HIType');
    optsStoc = rmfield(optsStoc,'stocName');
    optsStoc = rmfield(optsStoc,'noise');


    if(~fixedInitial)

        xpStruct.xEq = stocStruct(1).xInitial;
        xpStruct.pEq = stocStruct(1).pIF;

        opts.optsPlot.geom        = optsPlot.geom;
        opts.optsPlot.nBins       = optsPlot.nBins;
        opts.optsPlot.dim         = optsPlot.dim;
        opts.optsPlot.plotTimes   = optsPlot.plotTimes;
        opts.optsPlot.nParticlesS = optsPlot.nParticlesS;
        opts.optsPlot.mS          = optsPlot.mS;
        opts.optsPlot.rMin        = optsPlot.rMin;
        opts.optsPlot.rMax        = optsPlot.rMax;

        opts.optsStoc = optsStoc;
        opts.optsPhys = optsPhys;
        
%        redoPlotData = ~loadSamples;
          redoPlotData = true;

        fprintf(1,'Calculating Initial equilibria values ... ');
        eqDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Initial' filesep 'Equilibria'];
        initialEq = DataStorage(eqDir,@getEquilibria,opts,xpStruct,redoPlotData);
        fprintf(1,'Finished\n');

    else
        initialEq=[];
    end
    
    %--------------------------------------------------------------------------
    % get final equilibrium quantities to add to plots
    %--------------------------------------------------------------------------
    
    if(optsStocFull(1).sampleFinal)
    
        xpStruct.xEq = stocStruct(1).xFinal;
        xpStruct.pEq = stocStruct(1).pIF;

        opts.optsPlot.geom        = optsPlot.geom;
        opts.optsPlot.nBins       = optsPlot.nBins;
        opts.optsPlot.dim         = optsPlot.dim;
        opts.optsPlot.plotTimes   = optsPlot.plotTimes;
        opts.optsPlot.nParticlesS = optsPlot.nParticlesS;
        opts.optsPlot.mS          = optsPlot.mS;
        opts.optsPlot.rMin        = optsPlot.rMin;
        opts.optsPlot.rMax        = optsPlot.rMax;

        opts.optsStoc = optsStoc;
        opts.optsPhys = optsPhys;

%        redoPlotData = ~loadSamples;
         redoPlotData = true;
        
        fprintf(1,'Calculating Final equilibria values ... ');
        eqDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Final' filesep 'Equilibria'];
        finalEq = DataStorage(eqDir,@getEquilibria,opts,xpStruct,redoPlotData);
        fprintf(1,'Finished\n');

    else
        finalEq=[];
    end

    equilibria(1).data=initialEq;
    equilibria(2).data=finalEq;
    
    %--------------------------------------------------------------------------
    % get dynamic quantities (rho, v, means, etc)
    %--------------------------------------------------------------------------
    
    if(anyPlots) % ony if doing any non-particle plots
    
        opts.optsPhys = optsPhys;

        opts.optsPlot.geom        = optsPlot.geom;
        opts.optsPlot.nBins       = optsPlot.nBins;
        opts.optsPlot.dim         = optsPlot.dim;
        opts.optsPlot.plotTimes   = optsPlot.plotTimes;
        opts.optsPlot.nParticlesS = optsPlot.nParticlesS;
        opts.optsPlot.mS          = optsPlot.mS;
        if(optsPhys.dim == 2)
            opts.optsPlot.fixedBins   = optsPlot.fixedBins;
            opts.optsPlot.rMin        = optsPlot.rMin;
            opts.optsPlot.rMax        = optsPlot.rMax;
        end


        nStoc = length(optsStocFull);

        fprintf(1,'Calculating dynamics values ... ');
        for iStoc = 1:nStoc
            optsStoc = optsStocFull(iStoc);

            redoPlotData = ~optsStoc.loadStoc;
    %         redoPlotData = true;

            optsStoc = rmfield(optsStoc,'loadSamples');
            optsStoc = rmfield(optsStoc,'sampleFinal');
            optsStoc = rmfield(optsStoc,'poolsize');
            optsStoc = rmfield(optsStoc,'doStoc');
            optsStoc = rmfield(optsStoc,'loadStoc');
            optsStoc = rmfield(optsStoc,'saveStoc');
            optsStoc = rmfield(optsStoc,'stocName');

            opts.optsStoc = optsStoc;

            plotDataDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Dynamics' filesep 'Plotting'];
            stocPlotStruct(iStoc) = DataStorage(plotDataDir,@getAverages,opts,stocStruct(iStoc),redoPlotData);

        end
        fprintf(1,'Finished\n');
    
    end

else
    
    stocPlotStruct = [];
    equilibria(1).data = [];
    equilibria(2).data = [];
    
end

if(~isempty(DDFTStruct))
    
    nDDFT = length(DDFTStruct);
    
    for iDDFT = 1:nDDFT
                
        if(isa(DDFTStruct(1),'DDFT_2D'))

            rhoI = DDFTStruct(iDDFT).dynamicsResult.rho_t(:,:,1);
            rhoF = DDFTStruct(iDDFT).dynamicsResult.rho_t(:,:,end);

            DDFTStruct(iDDFT).dynamicsResult.rhoI = rhoI;
            DDFTStruct(iDDFT).dynamicsResult.rhoF = rhoF;

            DDFTAveragesStruct      = getAveragesDDFT([],DDFTStruct(iDDFT).dynamicsResult);
            DDFTStruct(iDDFT).dynamicsResult.rMean    = DDFTAveragesStruct.rMean;
            DDFTStruct(iDDFT).dynamicsResult.fluxMean = DDFTAveragesStruct.fluxMean;
            
        else

            rhoI = DDFTStruct(iDDFT).rho_t(:,:,1);
            rhoF = DDFTStruct(iDDFT).rho_t(:,:,end);

            DDFTStruct(iDDFT).rhoI = rhoI;
            DDFTStruct(iDDFT).rhoF = rhoF;

            DDFTAveragesStruct      = getAveragesDDFT([],DDFTStruct(iDDFT));
            DDFTStruct(iDDFT).rMean    = DDFTAveragesStruct.rMean;
            DDFTStruct(iDDFT).fluxMean = DDFTAveragesStruct.fluxMean;
            
        end
          
    end
    
    DDFTPlotStruct = DDFTStruct;
    
else
   
    DDFTPlotStruct = [];
    
end

if(anyPlots)
    plotDir = optsPlot.plotDir;
elseif(anyPlotsP)
    plotDir = optsPlotParticles.plotDir;
else
    plotDir = [];
end

if(anyPlots || anyPlotsP)
    if(~exist(plotDir,'dir'))
        mkdir(plotDir);
    end
end

if(anyPlots)

    %--------------------------------------------------------------------------
    % make movie
    %--------------------------------------------------------------------------

    if(optsPlot.doMovieGif || optsPlot.doMovieAvi || optsPlot.doMovieSwf || optsPlot.doPdfs) 
        movieFile = makeMovie(stocPlotStruct,DDFTPlotStruct,optsPlot,optsPhys,equilibria(2));
        plotFiles = cat(2,plotFiles,movieFile);
    end

    %--------------------------------------------------------------------------
    % make intial and final plots
    %--------------------------------------------------------------------------

    if(optsPlot.doInitialFinal)
        IFFile = plotInitialFinal(stocPlotStruct,DDFTPlotStruct,optsPlot,equilibria);
        plotFiles = cat(2,plotFiles,IFFile);
    end

    %--------------------------------------------------------------------------
    % make equilibrium plots
    %--------------------------------------------------------------------------

    if(optsPlot.doEquilibria)
        eqFile = plotEquilibria(stocStruct,DDFTStruct,optsPlot,optsPhys,equilibria);
        plotFiles = cat(2,plotFiles,eqFile);
    end

    %--------------------------------------------------------------------------
    % make mean plot
    %--------------------------------------------------------------------------

    if(optsPlot.doMeans)
        %pdfFile=fileStruct.plotFile{1};
        meanFile = plotMeans(stocPlotStruct,DDFTPlotStruct,optsPlot,equilibria);
        plotFiles = cat(2,plotFiles,meanFile);
    end

    %--------------------------------------------------------------------------
    % make error snapshot plots
    %--------------------------------------------------------------------------

    if(optsPlot.doSnapshotsError)
        errorFile = plotSnapshotsError2D(stocPlotStruct,DDFTPlotStruct,optsPlot,equilibria);
        plotFiles = cat(2,plotFiles,errorFile);
    end

    
end    
   
%--------------------------------------------------------------------------
% get mean (over runs) position and momentum of each particle at each time.
% only useful when particles are 'distinguishable'
%--------------------------------------------------------------------------

if(anyPlotsP)
    meanStruct=getMeansParticles(stocStruct);
end

% 
% %--------------------------------------------------------------------------
% % make intial and final particle plots
% %--------------------------------------------------------------------------
% 
% if(optsStruct.doInitialFinalP)
%     gifFile=fileStruct.plotFileP{2};
%     plotInitialFinalParticles(meanStruct,optsPlotParticles,gifFile);
% end
% 
% %--------------------------------------------------------------------------
% % make particle movie
% %--------------------------------------------------------------------------
% 

if(anyPlotsP)
    if(optsPlotParticles.doMovieGifP || optsPlotParticles.doMovieSwfP || optsPlotParticles.doPdfsP)
        makeMovieParticles(meanStruct,optsPlotParticles)
    end
end
    
%--------------------------------------------------------------------------
% make custom plot
%--------------------------------------------------------------------------

% if(optsStruct.doCustom)
%     path('CustomPlots',path)
%     customPlot=str2func(optsStruct.custom);
%     customPlot(stocStruct,DDFTStruct,equilibria,optsPlot,optsPhys);
%     rmpath('CustomPlots');
% end
% 
% %--------------------------------------------------------------------------
% % make custom particle plot
% %--------------------------------------------------------------------------
% 
% if(optsStruct.doCustomP)
%     path('CustomPlots',path)
%     customPlot=str2func(optsStruct.custom);
%     customPlot(meanStruct,optsPlotParticles);
%     rmpath('CustomPlots');
% end

end