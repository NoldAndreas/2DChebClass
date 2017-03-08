function outputFile = plotEquilibria1D(stoc,ddft,optsPlot,equilibria)
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


fprintf(1,'Making equilibria plots ... ');

% number of stochastic and DDFT plots to do
nStoc=size(stoc,2);
nDDFT=size(ddft,2);


% easy way to check if any of the computations include momentum
doP=~isempty(optsPlot.legTextP);

% set up figure
fullscreen = get(0,'ScreenSize');
hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
% set background colour to white
set(hRPf,'Color','w');

% if any of the calculations contain momentum, plot them to a visible axis,
% otherwise create a dummy set of axes which aren't visible
if(doP)    
    handles=tightsubplot(2,1,0.075,0.05,0.05);
    hRa=handles(1);
    hPa=handles(2);
else
    hRa=axes;
    % invisible axes off the figure
    hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
end

% figure and axes handles to be passed to plotting functions
handlesRP=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa);

% choose appropriate text for axis labels
if(optsPlot.symbolLabels)
    rhoText='$\rho \;\;$';
    %vText='$p \;\;$';     
    vText='Flux';     
    rText='$r$';
else
    geom=optsPlot.geom;
    if(strcmp(geom,'spherical'))
        %vText='Radial velocity';  
        vText='Radial flux';  
        rText='Radius';
    else
        %vText='Velocity';  
        vText='Flux';  
        rText='Distance';
    end
    rhoText='Number of particles';
end

%----------------------------------------------------------------------
% Get time, colours and types
%----------------------------------------------------------------------

% get time to add to plot
plotTimes=optsPlot.plotTimes;
plotTime=plotTimes(1);

% surface colours and types
if(nStoc>0)
     optsPlot.faceColour=optsPlot.eqColour{1};
     % get type info -- used to decide whether to plot the velocity
     stocType=optsPlot.stocType;
end

if(nDDFT>0)
    optsPlot.faceColour=optsPlot.eqColour{1};
    % get type info -- used to decide whether to plot the velocity
    DDFTType=optsPlot.DDFTType;
end

% and file to save in
outputFile = optsPlot.eqFile;


if(nStoc>0)

    rho   = equilibria(1).data.REq;
    v     = equilibria(1).data.vEq;
    boxes = equilibria(1).data.xEq;
    
    optsPlot.type=stocType(1,:);
    optsPlot.lineStyle = '-';

    optsPlot.lineColour = optsPlot.lineColourStoc{1};
    optsPlot.lineStyle = optsPlot.lineStyleStoc{1};
    optsPlot.lineMarker = optsPlot.lineMarkerStoc{1};

    
    plotRhoVdistStoc(rho,v,boxes,optsPlot,handlesRP,stocType(1));

    hold(hRa,'on');
    hold(hPa,'on');
 
end


if(nDDFT>0)

    % get rho and flux
    rho=ddft(1).rho_t;
    flux=ddft(1).flux_t;
   
    optsPlot.type=DDFTType(1,:);
    
    optsPlot.lineColour = optsPlot.lineColourDDFT{1};
    optsPlot.lineStyle = optsPlot.lineStyleDDFT{1};
    optsPlot.lineMarker = optsPlot.lineMarkerDDFT{1};

    % get values at appropriate time
    rhot  = rho(:,:,1);
    fluxt = flux(:,:,1);

    Interp = ddft(1).shape.Interp;
    
    plotRhoVdistDDFT(rhot,fluxt,Interp,optsPlot,handlesRP,optsPlot.type);

    hold(hRa,'on');
    hold(hPa,'on');

end

     

%----------------------------------------------------------------------
% Set axes, legend, add time
%----------------------------------------------------------------------

% and x ranges
rMin=optsPlot.rMin;
rMax=optsPlot.rMax;
pMin=optsPlot.pMin;
pMax=optsPlot.pMax;
% and y ranges
RMin=optsPlot.RMin;
RMax=optsPlot.RMax;
PMin=optsPlot.PMin;
PMax=optsPlot.PMax;
% get legend settings
legPos=optsPlot.legPos;
legTextR=optsPlot.legTextR;
legTextP=optsPlot.legTextP;

plotTime = [];

fixPlot(hRa,rMin,rMax,RMin,RMax,rText,rhoText,plotTime,legPos,legTextR);

if(~isempty(legTextP))  
    fixPlot(hPa,pMin,pMax,PMin,PMax,rText,vText,plotTime,legPos,legTextP);
end

% write the figure files
save2pdf(outputFile,hRPf);
%close(hRPf);

fprintf(1,'Finished\n');
