function outputFiles = plotSnapshotsDDFTWall2D(stoc,ddft,optsPlot,equilibria)
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
plotPos = [1;ceil(nSaves/2);ceil(nSaves)]; %away

yMin = 0;
yMax = 0.4;  % towards and away


nComps = nDDFT;

lineStyles = {':','-','--'};


for iPlot=1:nPlots  

    %----------------------------------------------------------------------
    % Set up figure
    %----------------------------------------------------------------------

    hf=figure('Position',[100,100,800,800]);
    ha = axes;
    set(hf,'Color','w');
    
    for iComp = 1:nComps
    
        
        %----------------------------------------------------------------------
        % DDFT data plots
        %----------------------------------------------------------------------

        
        IDC = ddft(iComp).IDC;

        PlotArea.y1Min = max(optsPlot.rMin(1),min(IDC.Pts.y1));
        PlotArea.y2Min = max(optsPlot.rMin(2),min(IDC.Pts.y2));
        PlotArea.y1Max = min(optsPlot.rMax(1),max(IDC.Pts.y1));
        PlotArea.y2Max = min(optsPlot.rMax(2),max(IDC.Pts.y2));
        PlotArea.N1 = 50;
        PlotArea.N2 = 50;

        
        InterpRho = IDC.InterpolationPlotCart(PlotArea,false);
        x1=InterpRho.pts1;
        x2=InterpRho.pts2;
        
        % get rho, v and r values
        rho=ddft(iComp).dynamicsResult.rho_t;

        rhot = rho(:,:,plotPos(iPlot));
        rhot = real(InterpRho.InterPol*rhot);
        
        x2Min = min(x2);
        x2MinPos = (x2==x2Min);
        
        rhoWall = rhot(x2MinPos);
        
        lineStyle = lineStyles{iComp};
        lineColour = optsPlot.lineColourDDFT{iComp};
        plot(ha,x1(x2MinPos),rhoWall,'Color',lineColour{1},'lineStyle',lineStyle);
        

        hold(ha,'on');
        


        %----------------------------------------------------------------------
        % Set axes, legend, add time
        %----------------------------------------------------------------------


        optsPlot.xLab='x';
        optsPlot.yLab='Density at wall';

    end % for iComp
    
    ylim([yMin,yMax]);
    
    xLims = get(ha,'XLim');
    yLims = get(ha,'YLim');
    
    textX=xLims(1) + 0.8*(xLims(2)-xLims(1));
    textY=yLims(1)+0.9*(yLims(2)-yLims(1));
    
    plotTimes=optsPlot.plotTimes;
    plotTime=plotTimes(plotPos(iPlot));
    
    text(textX,textY,['t=' num2str(plotTime, '%6.4f')]);

    set(get(ha,'XLabel'),'String','$x$','Interpreter','LaTex','FontSize',20)
    set(get(ha,'YLabel'),'String','Density at wall','Interpreter','LaTex','FontSize',20)

    outputFile = [optsPlot.plotDir filesep 'snapshotWall-' num2str(plotTime,3) '.fig'];
                            
    savefig(hf,outputFile);
    
end % for iPlot

fprintf(1,'Finished\n');
