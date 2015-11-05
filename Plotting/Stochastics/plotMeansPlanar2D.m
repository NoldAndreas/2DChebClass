function outputFile = plotMeansPlanar2D(stoc,ddft,optsPlot,equilibria)
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

% geometry and dimension needed to calculate means
geom=optsPlot.geom;
dim=optsPlot.dim;

% pre-calculate means for both stochastic and DDFT data
% meanrStoc=zeros(nPlots,nSpecies,dim,nStoc);
% meanvStoc=meanrStoc;
% 
% for iStoc=1:nStoc
%     % loop through each set of stochastic data
%     x=stoc(iStoc).x;
%     p=stoc(iStoc).p;
%     
%     % compute means for each time, species and stochastic calculation
%     [meanrStoc(:,:,:,iStoc),meanvStoc(:,:,:,iStoc)]=getRVmeansStoc2D(x,p,geom,nParticlesS,mS,optsPlot.saveFileStocMeans{iStoc});
% end

if(nStoc>0)
    meanrStoc = stoc.rMean;
    meanvStoc = stoc.vMean;
end

meanrDDFT=zeros(nPlots,nSpecies,dim,nDDFT);
meanvDDFT=meanrDDFT;

for iDDFT=1:nDDFT 
    % loop through each set of DDFT data and compute / read in means
    [meanrDDFT(:,:,:,iDDFT),meanvDDFT(:,:,:,iDDFT)]=getRVmeansDDFT2D(ddft(iDDFT),geom);
end

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
hMRa=cell(1,dim);
hMPa=cell(1,dim);

if(~isempty(legTextP))    
    %handles=tightsubplot(2,1,0.075,0.075,0.075);
    handles=tightsubplot(2,dim,0.025,0.075,0.075);
    for iDim=1:dim
        hMRa{iDim}=handles(iDim);
        hMPa{iDim}=handles(dim+iDim);
    end
else
    handles=tightsubplot(1,dim,0.075,0.075,0.075);
    hDummy=axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
    for iDim=1:dim
        hMRa{iDim}=handles(iDim);
        hMPa{iDim}=hDummy;
    end
end

% figure and axis handles to be passed to plotting functions
handlesM=struct('hRPf',hRPf,'hRa',hMRa,'hPa',hMPa);

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
    
    meanrStocS=meanrStoc(:,:,:,iStoc);
    meanvStocS=meanvStoc(:,:,:,iStoc);
    
    % plot means of r and v up to current time    
    for iDim=1:dim
        plotRVmeans(meanrStocS(:,:,iDim),meanvStocS(:,:,iDim),plotTimes,nPlots,optsPlot,handlesM(iDim),stocType(iStoc));
        hold(hMRa{iDim},'on');
        hold(hMPa{iDim},'on');
    end
end

%----------------------------------------------------------------------
% DDFT plotting
%----------------------------------------------------------------------

for iDDFT=1:nDDFT
    % choose line style
    optsPlot.lineStyle=lineStyleDDFT{iDDFT};
    optsPlot.lineMarker=lineMarkerDDFT{iDDFT};
    optsPlot.lineColour=lineColourDDFT{iDDFT};

    meanrDDFTS=meanrDDFT(:,:,:,iDDFT);
    meanvDDFTS=meanvDDFT(:,:,:,iDDFT);
    
    for iDim=1:dim
        % plot means of r and v up to current time
        plotRVmeans(meanrDDFTS(:,:,iDim),meanvDDFTS(:,:,iDim),plotTimes,nPlots,optsPlot,handlesM(iDim),DDFTType(iDDFT));
        hold(hMRa{iDim},'on');
        hold(hMPa{iDim},'on');
    end
end

%----------------------------------------------------------------------
% Sort out axes and legend
%----------------------------------------------------------------------

% y axis limits for position plot
RMMin=optsPlot.RMMin;
RMMax=optsPlot.RMMax;

switch optsPlot.geom
    case 'planar2D'
        yLabelsR={'$\bar x  \;\;$','$\bar y  \;\;$'};
        yLabelsP={'$\bar v_x  \;\;$','$\bar v_y  \;\;$'};
    case 'polar2D'
        yLabelsR={'$\bar r  \;\;$','$\bar \theta  \;\;$'};
        yLabelsP={'$\bar v_r  \;\;$','$\bar v_\theta  \;\;$'};
end

% modify position plot
for iDim=1:dim
    %fixPlot(hMRa{iDim},tMin,tMax,RMMin,RMMax,'$t$','$\bar r  \;\;$',[],legPos,legTextR)    
    fixPlot(hMRa{iDim},tMin,tMax,RMMin(iDim),RMMax(iDim),'$t$',yLabelsR{iDim},[],legPos,legTextR)
end

% only modify momentum axis if we've plotted anything there
if(~isempty(legTextP))
    % y axis limits for momentum plot
    PMMin=optsPlot.PMMin;
    PMMax=optsPlot.PMMax;
    %modify momentum plot
    for iDim=1:dim
        %fixPlot(hMPa{iDim},tMin,tMax,PMMin,PMMax,'$t$','$\bar p  \;\;$',[],legPos,legTextP)
        fixPlot(hMPa{iDim},tMin,tMax,PMMin(iDim),PMMax(iDim),'$t$',yLabelsP{iDim},[],legPos,legTextP)
    end
    % remove tick marks and label from position plot
    set(hMRa,'XTickLabel',[]);
    %set(get(hMRa,'XLabel'),'String','');
end

% add equilibrium dotted lines
nEq=size(equilibria,2);
for iEq=1:nEq
    if(~isempty(equilibria(iEq).data))
        for iDim=1:dim
            data.rMeanEq=equilibria(iEq).data.rMeanEq(:,:,iDim);
            data.vMeanEq=equilibria(iEq).data.vMeanEq(:,:,iDim);
            data.plotTimes=equilibria(iEq).data.plotTimes;
            %addEqLines([],handlesM,equilibria(iEq).data,optsPlot);
            addEqLines([],handlesM(iDim),data,optsPlot);
        end
    end
end

if(nDDFT>0)
%   addEqLinesDDFT([],handlesM,ddft,optsPlot,true,true);
end

outputFile = optsPlot.meansFile;

save2pdf(outputFile,hRPf);

fprintf(1,'Finished\n');