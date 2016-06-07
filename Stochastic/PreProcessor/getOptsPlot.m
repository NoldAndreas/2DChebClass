function [optsPlot,optsPlotParticles]=getOptsPlot(optsStruct)
% [optsPlotGIF,optsPlotParticles]=getOptsPlot(optsStruct,fileStruct)
%   Collates plotting options for standard plots and particle plots
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%    FOR getOptsPlot
%         plotTimes         save times (nPlots,1)
%         geom              geometry (string)
%         stocDim           stochastic (full) dimension (integer)
%         nParticlesS       number of particles (nSpecies,1)
%         mS                mass of particles (nSpecies,1)
%                     
%         rMin              position axes min {dim,1}
%         rMax              position axes max {dim,1}
%         pMin              velocity axes min {dim,1} 
%         pMax              velocity axes max {dim,1} 
%         RMin              rho plot min, scalar or {rG,rI,rF} 
%         RMax              rho plot max, scalar or {rG,rI,rF}            
%         PMin              v plot min, scalar or {vG,vI,vF} 
%         PMax              v plot max, scalar or {vG,vI,vF}            
%         RMMin             mean position min {dim,1}
%         RMMax             mean position max {dim,1}                   
%         PMMin             mean momentume min {dim,1}
%         PMMax             mean momentume min {dim,1}
%                     
%         nBins             number of bins for histogramming (dim,1)    
%         fixedBins         t/f whether to use fixed bins (or adaptive)
%
%         plotDensity       t/f plot density rather than distribution
%                            (used for spherical geom)                    
%         plotCurrent       t/f plot current (rho v) rather than v                   
%         plotType          'surf' or 'contour' for 2D plots
%         plotEquilibria    t/f add equilibria lines  
%         vCutoff           scalar doesn't plot v for areas with rho< vCutoff*max(rho)                  
%         viewPoint         [Az,El] for view angle of 3D plots                   
%         symbolLabels      t/f whether to use symbols (t) or words (f) for labels
% 
%         doMovieGif        t/f do gif movie
%         doMovieSwf        t/f do flash movie
%         doPdfs            t/f save pdfs for movie frames
% 
%         dpi               dpi for rendering
%         fps               fps for movies
%         bitmap            t/f bitmap render in swf
%         quiet             t/f supress gs output
%
%         anyPlots          t/f do any standard plots
%         anyPlotsP         t/f do any particle plots
%         anyDDFT           t/f do any DDFT calculations
%
%    FOR makeLegLine
%         doStoc       whether to do each of the stochastic simulations, 
%                        of the form {true, false,...}
%         doDDFT       same for DDFT
%         sampleFinal  true/false whether we sample, and hence plot, 
%                        the final eq
%         stocName     legend labels for each stochastic simulation, 
%                        of the form {name1,name2,...}
%         DDFTName     same for DDFT
%         stocStyle    line styles for plotting stochastic data, 
%                        of the same form
%         DDFTStyle    same for DDFT
%         stocColour   line colours for plotting stochastic data, 
%                        of the same form
%         DDFTColour   same for DDFT
%         stocType     type of stochastic simulations, {'rv','r',...}
%         DDFTType     same for DDFT
%         legPos       legend position, e.g. 'NorthWest'
%         oneLeg       either 'off' or anything else, to plot a single 
%                        legend at the top of the figure
%         perRow       integer max no. labels per row in single legend
%       
%   fileStruct  -- struct of size [1 1] containing
%                   pdfDir:        string save directory for pdf files for
%                                   movies
%                   movieFile:     string save file for gif/swf movies
%
%  PARTICLE PLOT OPTIONS CURRENTLY UNDOCUMENTED!!
%
% OUTPUTS:
%   optsPlotGIF        -- struct of size [1 3] containing all the above data
%                          required for standard plots
%   optsPlotParticles  -- struct of size [1 3] containing all the above data
%                          required for particle plots

%--------------------------------------------------------------------------
% Get legend data
%--------------------------------------------------------------------------

optsLeg=makeLegLine(optsStruct);

%--------------------------------------------------------------------------
% Get potential plotting options
%--------------------------------------------------------------------------

%optsPlotV=getOptsPlotV(optsStruct);

%--------------------------------------------------------------------------
% Construct standard plotting options
%--------------------------------------------------------------------------

global dirData

plotDir    = [dirData filesep optsStruct.potNames filesep 'Output'];
pdfDir     = [plotDir filesep 'pdfs' filesep];
movieFile  = [plotDir filesep 'dynamics'];
meansFile  = [plotDir filesep 'means.pdf'];
eqFile     = [plotDir filesep 'equlibria.pdf'];
IFFiles    = {{[plotDir filesep 'initial.pdf'], [plotDir filesep 'final.pdf']}};

if(optsStruct.anyPlots || optsStruct.anyDDFT)  

    optsPlot=struct('plotTimes',optsStruct.plotTimes, ...
                    'geom',optsStruct.geom, ...
                    'dim',optsStruct.stocDim, ...    
                    'nParticlesS',optsStruct.nParticlesS, ...
                    'mS',optsStruct.mS, ...
                    'rMin',optsStruct.rMin, 'rMax',optsStruct.rMax, ...
                    'pMin',optsStruct.pMin, 'pMax',optsStruct.pMax, ...
                    'RMin',optsStruct.RMin, 'RMax',optsStruct.RMax, ...
                    'PMin',optsStruct.PMin, 'PMax',optsStruct.PMax, ...
                    'RMMin',optsStruct.RMMin,'RMMax',optsStruct.RMMax, ...
                    'PMMin',optsStruct.PMMin,'PMMax',optsStruct.PMMax, ...
                    'plotDensity',optsStruct.plotDensity, ...
                    'plotCurrent',optsStruct.plotCurrent, ...
                    'plotEquilibria',optsStruct.plotEquilibria, ...
                    'viewPoint',optsStruct.viewPoint, ...
                    'doMovieGif',optsStruct.doMovieGif, ...
                    'doMovieAvi',optsStruct.doMovieAvi, ...
                    'doMovieSwf',optsStruct.doMovieSwf, ...
                    'doPdfs',optsStruct.doPdfs, ...
                    'doInitialFinal',optsStruct.doInitialFinal,...
                    'doEquilibria',optsStruct.doEquilibria,...
                    'doMeans',optsStruct.doMeans,...
                    'doSnapshotsError',optsStruct.doSnapshotsError,...
                    'doSnapshotsDDFT',optsStruct.doSnapshotsDDFT,...
                    'dpi',optsStruct.dpi, ...
                    'fps',optsStruct.fps, ...
                    'bitmap',optsStruct.bitmap, ...
                    'quiet',optsStruct.quiet, ...
                    'symbolLabels',optsStruct.symbolLabels, ...
                    'plotType',optsStruct.plotType, ...
                    'separateSpecies',optsStruct.separateSpecies, ...
                    'separateError',optsStruct.separateError, ...
                    'separateComp',optsStruct.separateComp, ...
                    'plotDir',plotDir,'pdfDir',pdfDir, 'movieFile',movieFile, ...
                    'meansFile',meansFile,'IFFiles',IFFiles,'eqFile',eqFile);

    % merge in the legend and line options            
    %optsPlot=mergeStruct(optsPlot,optsLeg,optsPlotV,optsStruct.potParamsPlot);
    optsPlot=mergeStruct(optsPlot,optsLeg,optsStruct.potParamsPlot);
    
    if(optsStruct.anyStoc)
        optsPlotStocGIF=struct('nBins',optsStruct.nBins, ...
                    'fixedBins',optsStruct.fixedBins, ...
                    'vCutoff',optsStruct.vCutoff );
        optsPlot=mergeStruct(optsPlot,optsPlotStocGIF);
    end
    
else 
    optsPlot=[];
end

%--------------------------------------------------------------------------
% Construct particle plotting options
%--------------------------------------------------------------------------

if(optsStruct.anyPlotsP)        
    optsPlotParticles=struct('shift',optsStruct.shift, ...
                'viewPoint',optsStruct.viewPoint, ...
                'followType',optsStruct.followType, ...
                'comPlot',optsStruct.comPlot, ...
                'renormalize',optsStruct.renormalize, ...
                'relabel',optsStruct.relabel, ...
                'lims',optsStruct.lims, ...
                'ticks',optsStruct.ticks, ...
                'labs',optsStruct.labs, ...
                'colours',optsStruct.colours, ...
                'bath',optsStruct.bath, ...
                'sigma',optsStruct.sigma, ...
                'pScale',optsStruct.pScale, ...
                'plotTimes',optsStruct.plotTimes, ...
                'doMovieGifP',optsStruct.doMovieGifP, ...
                'doMovieSwfP',optsStruct.doMovieSwfP, ...
                'doPdfsP',optsStruct.doPdfsP, ...
                'dpi',optsStruct.dpi, ...
                'fps',optsStruct.fps, ...
                'bitmap',optsStruct.bitmap, ...
                'quiet',optsStruct.quiet, ...
                'plotDir',plotDir,'pdfDir',pdfDir, 'movieFile',movieFile);
              
else
    optsPlotParticles=[];
end

end