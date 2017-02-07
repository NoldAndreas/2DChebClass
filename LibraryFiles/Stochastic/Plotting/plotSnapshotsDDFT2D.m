function outputFiles = plotSnapshotsDDFT2D(stoc,ddft,optsPlot,equilibria)
%plotInitialFinal(stoc,ddft,optsPlotGIF,xInitial,xFinal,pEq,pdfFile)
%   makes initial and final plots from given stochastic and DDFT data
%
% INPUTS: 
%  stoc -- a structure of size [1 nStoc] of the form (for nStoc=2)
%           stoc=struct('x',{x1,x2},'p',{p1,p2}) with
%             x   ( (nRuns,dim*nParticles,nTimes) matrix of positions)
%             p   ( (nRuns,dim*nParticles,nTimes) matrix of momenta)
%  ddft -- a structure of size [1 nDDFT] of the form (for nDDFT=2)
%           ddft=struct('rhov',{rhov1,rhov2},'r',{r1,r2},'w',{w1,w2}):
%             rhov  ( (2*N,nTimes) matrix with columns [rho; v] specified at 
%                   grid points r and each of the nTimes)
%             r     ( (N,1) vector of grid points)
%             w     ( (1,N) vector of corresponding weights)
%
%  optsPlotGIF -- a structure of size [1 3] (General, Initial, Final) 
%                  containing:
%             plotTimes       (vector of plot times)
%             rMin            (minimum x value for density plots)
%             rMax            (maximum x value for density plots)
%             pMin            (minimum x value for momentum plots)
%             pMax            (maximum x value for momentum plots)
%             RMin            (minimum y value for density plots)
%             RMax            (maximum y value for density plots)
%             PMin            (minimum y value for momentum plots)
%             PMax            (maximum y value for momentum plots)
%             RMMin           (minimum y value for mean position plots)
%             RMMax           (maximum y value for mean position plots)
%             PMMin           (minimum y value for mean momentum plots)
%             PMMax           (maximum y value for mean momentum plots)
%             lineStyleStoc   (linestyles for stochastic plots, should be 
%                               of the form ['b  ';'--b'])
%             lineStyleDDFT   (linestyles for DDFT plots, should be
%                               of the form ['b  ';'--b'])
%             legPos          (legend position, e.g. 'SouthEast')
%             oneLeg          (not equal to 'off' to just plot one legend)
%             perRow          (number of lines per row in oneLeg)
%             legText         (legend text, should be of the form 
%                               {'label1' 'label2'})
%             geom            (geometry; either 'planar' or 'spherical')
%             dim             (dimension)
%
%  xInitial  --  (nSamples,dim*nParticles) matrix of initial positions
%  xFinal    --  (nSamples,dim*nParticles) matrix of final positions
%               OR [] for no final sampling
%  pEq       --  (nSamples,dim*nParticles) matrix of initial and final 
%                momenta
%
%  pdfFile    --     1x3 cell of strings for output files (only last two
%                     are used)


fprintf(1,'Making snapshot DDFT plots ... ');

% number of stochastic and DDFT plots to do
nDDFT=size(ddft,2);

% easy way to check if any of the computations include momentum
%doP=~isempty(optsPlot.legTextP);

outputFiles = {};

nPlots = 3;
nSaves = length(optsPlot.plotTimes);

% this cuts down the final plotting time!!
%nSaves = ceil(nSaves/2);
%nSaves = ceil(nSaves/4);

%dPlot = floor((nSaves - 1)/(nPlots-1));
%plotPos = 1:dPlot:nSaves;


%plotPos = [1;ceil(nSaves/3);ceil(2*nSaves/3)]; % wall towards
plotPos = [1;ceil(nSaves/2);ceil(nSaves)]; %unbounded, away

nComps = nDDFT;

for iPlot=1:nPlots  

    for iComp = 1:nComps
    
        %----------------------------------------------------------------------
        % Set up figure
        %----------------------------------------------------------------------

        hSf=figure('Position',[100,100,800,800]);
        hSa = axes;
        hPSa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
        % set background colour to white
        set(hSf,'Color','w');

        hCf=figure('Position',[100,100,800,800]);
        hPCa = axes;
        hCa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
        % set background colour to white
        set(hCf,'Color','w');


        % figure and axes handles to be passed to plotting functions
        handlesS=struct('hRPf',hSf,'hRa',hSa,'hPa',hPSa);
        handlesC=struct('hRPf',hCf,'hRa',hCa,'hPa',hPCa); 

        %----------------------------------------------------------------------
        % Get time, colours and types
        %----------------------------------------------------------------------

        % get time to add to plot
        plotTimes=optsPlot.plotTimes;
        plotTime=plotTimes(plotPos(iPlot));
        
        lineColourDDFT=optsPlot.lineColourDDFT;
        lineStyleDDFT=optsPlot.lineStyleDDFT;
        % get type info -- used to decide whether to plot the velocity
        DDFTType=optsPlot.DDFTType;

        %----------------------------------------------------------------------
        % DDFT data plots
        %----------------------------------------------------------------------


        % get rho, v and r values
        rho=ddft(iComp).dynamicsResult.rho_t;

        optsPlot.type=DDFTType{iComp};

        % get values at appropriate time
        rhot=rho(:,:,plotPos(iPlot));

        if(strcmp(optsPlot.type,'rv'))
            v=ddft(iComp).v_IP;
            vt=v(:,:,plotPos(iPlot));
        else
            v=ddft(iComp).dynamicsResult.flux_t;
            fluxNorm = 0.1*max(max(max(v)));
            optsPlot.fluxNorm=fluxNorm;

            vt=v(:,:,plotPos(iPlot));
        end

        optsPlot.faceColour=lineColourDDFT{iComp};
        optsPlot.lineColour=lineColourDDFT{iComp};
        optsPlot.lineStyle=lineStyleDDFT{iComp};

        % plot the distributions
        plotRhoVdistDDFT2D(rhot,vt,ddft(iComp).IDC,optsPlot,handlesS);
        plotRhoVdistDDFT2D(rhot,vt,ddft(iComp).IDC,optsPlot,handlesC);

        hold(hSa,'on');
        hold(hPSa,'on');
        hold(hCa,'on');
        hold(hPCa,'on');

        


        %----------------------------------------------------------------------
        % Set axes, legend, add time
        %----------------------------------------------------------------------

        optsPlot.time=plotTime;

        optsPlot.xLab='x';
        optsPlot.yLab='z';

        optsPlot.xMin=optsPlot.rMin(1);
        optsPlot.xMax=optsPlot.rMax(1);
        optsPlot.yMin=optsPlot.rMin(2);
        optsPlot.yMax=optsPlot.rMax(2);

        optsPlot.zMin=optsPlot.RMin;
        optsPlot.zMax=optsPlot.RMax;

        optsPlot.zLab='Density';

        fixPlot2Dsurf(hSa,optsPlot);
        
        
        optsPlot.time = [];
        fixPlot2Dcontour(hPCa,optsPlot);

        outputFileS = [optsPlot.plotDir filesep 'snapshotSurface-' num2str(iComp) ...
                        '-' num2str(plotTime,3) '.fig'];
                    
        outputFileC = [optsPlot.plotDir filesep 'snapshotContour-' num2str(iComp) ...
                        '-' num2str(plotTime,3) '.fig'];
        
        savefig(hSf,outputFileS);
        savefig(hCf,outputFileC);
        
        % write the figure files
    %     save2pdf(outputFile,hRPf,100,true);
        %close(hRPf);

    %    outputFiles = cat(2,outputFiles,outputFile);

    end % for iComp
    
end % for iPlot

fprintf(1,'Finished\n');
