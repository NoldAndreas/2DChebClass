function plotMeans1D(stoc,ddft,optsPlot,equilibria)
% plotMeans(stoc,ddft,optsPlotGIF,equilibria,pdfFile)
%   plots mean position and momentum over time from given stochastic and DDFT data
%
% INPUTS: 
%  stoc -- a structure of size [1 nStoc] of the form (for nStoc=2)
%           stoc=struct('x',{x1,x2},'p',{p1,p2}) with
%             x   ( (nSamples,dim*nParticles,nTimes) matrix of positions)
%             p   ( (nSamples,dim*nParticles,nTimes) matrix of momenta)
%  ddft -- a structure of size [1 nDDFT] of the form (for nDDFT=2)
%           ddft=struct('rhov',{rhov1,rhov2},'r',{r1,r2},'w',{w1,w2}):
%             rhov  ( (2*N,nTimes) matrix with columns [rho; v] specified at 
%                   grid points r and each of the nTimes)
%             r     ( (N,1) vector of grid points)
%             w     ( (1,N) vector of corresponding weights)
%
%  optsPlot -- a structure of size [1 1] containing at least:
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
%                               ['b  ';'--b'])
%             lineStyleDDFT   (linestyles for DDFT plots, should be of the form
%                               of the form ['b  ';'--b'])
%             legPos          (legend position, e.g. 'SouthEast')
%             oneLeg          (not equal to 'off' to just plot one legend)
%             perRow          (number of lines per row in oneLeg)
%             legText         (legend text, should be of the form 
%                               {'label1' 'label2'})
%             geom            (geometry; either 'planar' or 'spherical')
%             dim             (dimension; either 1 or 3)
%
%  equilbria -- structure of size [1 nEq] containing
%             REq             (vector of rho at positions xR)
%             vEq             (vector of v at positions xR)
%             xR              (vector positions, centres of histogram bins)
%             rMeanEq         (vector of mean positions, 
%                               same length as plotTimes)
%             pMeanEq         (vector of mean momenta, 
%                               same length as plotTimes)
%             plotTimes       (vector of plot times)
%
%  pdfFile   -- string for output file name

fprintf(1,'Making mean plot ... ');

% number of stochastic and DDFT plots to do
nStoc=size(stoc,2);
nDDFT=size(ddft,2);

% number of species
nParticlesS=optsPlot.nParticlesS;
nSpecies=length(nParticlesS);

% species masses
mS=optsPlot.mS;

% get times
plotTimes=optsPlot.plotTimes;

% get initial and final times
tMin=plotTimes(1);
tMax=plotTimes(end);

% line styles and types
if(nStoc>0)
    lineStyleStoc=optsPlot.lineStyleStoc;
    lineMarkerStoc=optsPlot.lineMarkerStoc;
    lineColourStoc=optsPlot.lineColourStoc;
    % get type info -- used to decide whether to plot the velocity
    stocType=optsPlot.stocType;
else
    % choose line styles
    lineStyleStoc=[];
    lineMarkerStoc=[];
    lineColourStoc=[];
    % get type info -- used to decide whether to plot the velocity
    stocType=[];
end

if(nDDFT>0)
    lineStyleDDFT=optsPlot.lineStyleDDFT;
    lineMarkerDDFT=optsPlot.lineMarkerDDFT;
    lineColourDDFT=optsPlot.lineColourDDFT;
    % get type info -- used to decide whether to plot the velocity
    DDFTType=optsPlot.DDFTType;
else
    lineStyleDDFT=[];
    lineMarkerDDFT=[];
    lineColourDDFT=[];
    % get type info -- used to decide whether to plot the velocity
    DDFTType=[];
end

% get legend settings
legPos=optsPlot.legPos;
legTextR=optsPlot.legTextR;
legTextP=optsPlot.legTextP;

% number of plots in the movie
nPlots=length(plotTimes);

% set up figure
fullscreen = get(0,'ScreenSize');
if(fullscreen(3)>1300)
    fullscreen(3)=fullscreen(3)/2;
end

hRPf=figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
% set background colour to white
set(hRPf,'Color','w');

% if any of the calculations contain momentum, plot them to a visible axis,
% otherwise create a dummy set of axes which aren't visible
if(~isempty(legTextP))    
    %handles=tightsubplot(2,1,0.075,0.075,0.075);
    handles=tightsubplot(2,1,0.025,0.075,0.075);
    hMRa=handles(1);
    hMPa=handles(2);
else
    hMRa=axes;
    % invisible axes off the figure
    hMPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
end
    
% figure and axis handles to be passed to plotting functions
handlesM=struct('hRPf',hRPf,'hRa',hMRa','hPa',hMPa);

% if we're adding one legend, do it now
if(~strcmp(optsPlot.oneLeg,'off'))
    % combine legend labels, styles and colours, with the assumption
    % that they are all contained in legTextR (i.e. all are plotted in
    % the position plot
    optsPlot.legText=optsPlot.legTextR;
    optsPlot.lineStyle=cat(2,lineStyleStoc,lineStyleDDFT);
    optsPlot.lineMarker=cat(2,lineMarkerStoc,lineMarkerDDFT);
    optsPlot.lineColour=cat(2,lineColourStoc,lineColourDDFT);

    % create the legend
    addOneLeg(optsPlot,hRPf);
end 
   
%----------------------------------------------------------------------
% Stochastic plotting
%----------------------------------------------------------------------

for iStoc=1:nStoc
    % choose line style
    optsPlot.lineStyle=lineStyleStoc{iStoc};
    optsPlot.lineMarker=lineMarkerStoc{iStoc};
    optsPlot.lineColour=lineColourStoc{iStoc};
    
    plotRVmeans(stoc(iStoc).rMean,stoc(iStoc).vMean,plotTimes,nPlots,optsPlot,handlesM,stocType(iStoc));
    hold(hMRa,'on');
    hold(hMPa,'on');
end

%----------------------------------------------------------------------
% DDFT plotting
%----------------------------------------------------------------------

for iDDFT=1:nDDFT
    % choose line style
    optsPlot.lineStyle=lineStyleDDFT{iDDFT};
    optsPlot.lineMarker=lineMarkerDDFT{iDDFT};
    optsPlot.lineColour=lineColourDDFT{iDDFT};

    plotRVmeans(ddft(iDDFT).rMean,ddft(iDDFT).fluxMean,plotTimes,nPlots,optsPlot,handlesM,DDFTType(iDDFT));
    %plotRVmeans(meanrDDFTS,meanvDDFTS,plotTimes,nPlots,optsPlot,handlesM,DDFTType(iDDFT));
    hold(hMRa,'on');
    hold(hMPa,'on');

end

%----------------------------------------------------------------------
% Sort out axes and legend
%----------------------------------------------------------------------

% y axis limits for position plot
RMMin=optsPlot.RMMin;
RMMax=optsPlot.RMMax;
% modify position plot
fixPlot(hMRa,tMin,tMax,RMMin,RMMax,'$t$','$\bar r  \;\;$',[],legPos,legTextR)

% only modify momentum axis if we've plotted anything there
if(~isempty(legTextP))
    % y axis limits for momentum plot
    PMMin=optsPlot.PMMin;
    PMMax=optsPlot.PMMax;
    %modify momentum plot
    %fixPlot(hMPa,tMin,tMax,PMMin,PMMax,'$t$','$\bar p  \;\;$',[],legPos,legTextP)
    fixPlot(hMPa,tMin,tMax,PMMin,PMMax,'$t$','$\bar j  \;\;$',[],legPos,legTextP)
    
    % remove tick marks and label from position plot
    set(hMRa,'XTickLabel',[]);
    %set(get(hMRa,'XLabel'),'String','');
end

% add equilibrium dotted lines
nEq=size(equilibria,2);
for iEq=1:nEq
    if(~isempty(equilibria(iEq).data))
        addEqLines([],handlesM,equilibria(iEq).data,optsPlot)
    end
end

if(nDDFT>0)
   addEqLinesDDFT([],handlesM,ddft,optsPlot,true,true);
end

outputFile = optsPlot.meansFile;

save2pdf(outputFile,hRPf);

fprintf(1,'Finished\n');