function plotInitialFinal1D(stoc,ddft,optsPlot,equilibria)
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


fprintf(1,'Making initial and final plots ... ');

% number of stochastic and DDFT plots to do
nStoc=size(stoc,2);
nDDFT=size(ddft,2);

% if there are any stochastic plots, add in the initial and final data

%doEqStoc=[false;false;false];

if(nStoc>0)
    doEqStoc=[true;true];
else
    doEqStoc=[false;false];
end

if(nDDFT>0)
    doEqDDFT=[true;true];
else
    doEqDDFT=[false;false];
end

for iEq = 1:2
    if( ~isempty( equilibria(iEq).data ))
        doEqStoc(iEq)=true;
    end
    
    if( nDDFT>0 )
        doEqDDFT(iEq)=true;
    end
    
end

doP=~isempty(optsPlot.legTextP);


% choose appropriate text for axis labels
if(optsPlot.symbolLabels)
    rhoText='$\rho \;\;$';
    vText='$p \;\;$';     
    rText='$r$';
else
    geom=optsPlot.geom;
    if(strcmp(geom,'spherical'))
        vText='Radial velocity';  
        rText='Radius';
    else
        vText='Velocity';  
        rText='Distance';
    end
    rhoText='Number of particles';
end


% set initial and final time positions
plotPos(1)=1;
plotPos(2)=length(optsPlot.plotTimes);

for iPlot=1:2

    % get time to add to plot
    plotTimes=optsPlot.plotTimes;
    plotTime=plotTimes(plotPos(iPlot));

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
    
    % and file to save in
    outputFile=optsPlot.IFFiles{iPlot};
    
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
    
    % if we're adding one legend, do it now
    if(~strcmp(optsPlot.oneLeg,'off'))
        % create the legend
        addOneLeg(optsPlot,hRPf);
    end 
    
    %----------------------------------------------------------------------
    % Stochastic data plots
    %----------------------------------------------------------------------
        
    for iStoc=1:nStoc
    
        optsPlot.lineColour=lineColourStoc{iStoc};
        optsPlot.lineStyle=lineStyleStoc{iStoc};
        optsPlot.lineMarker=lineMarkerStoc{iStoc};
        
        rho   = stoc(iStoc).rho(:,:,plotPos(iPlot));
        v     = stoc(iStoc).v(:,:,plotPos(iPlot));
        boxes =  stoc(iStoc).boxes(:,:,plotPos(iPlot));
        
        optsPlot.plotTime=plotTime;
        
        plotRhoVdistStoc(rho,v,boxes,optsPlot,handlesRP,stocType(iStoc));
        %plotRhoVdistStoc(xt,pt,optsPlot,handlesRP,stocType(iStoc));
        
    end
    
    %----------------------------------------------------------------------
    % Equilibrium stochastic data plots
    %----------------------------------------------------------------------
    
    if(doEqStoc(iPlot))
        
%        addEqLines(handlesRP,[],equilibria(iPlot).data,optsPlot(iPlot))

    end
    
    %----------------------------------------------------------------------
    % DDFT data plots
    %----------------------------------------------------------------------
    
    for iDDFT=1:nDDFT
                
        % choose line style
        optsPlot.lineStyle=lineStyleDDFT{iDDFT};
        optsPlot.lineMarker=lineMarkerDDFT{iDDFT};
        optsPlot.lineColour=lineColourDDFT{iDDFT};
        
        % get rho and v
        rho=ddft(iDDFT).rho_t;
 %       v=ddft(iDDFT).v_t;
        
        % select correct time
        rho=rho(:,:,plotPos(iPlot));
%        v=v(:,:,plotPos(iPlot));
        v = zeros(size(rho));
        
        Interp = ddft(iDDFT).shape.Interp;
        
        % plot the distributions
        plotRhoVdistDDFT(rho,v,Interp,optsPlot,handlesRP,DDFTType(iDDFT));

        hold(hRa,'on');
        hold(hPa,'on');

    end
    
    if(doEqDDFT(iPlot))

        if(iPlot==2)
            doI=true;
            doF=false;
        else
            doI=false;
            doF=true;
        end
        
       % addEqLinesDDFT(handlesRP,[],ddft,optsPlot(iPlot),doI,doF)

    end
    
    %----------------------------------------------------------------------
    % Set axes, legend, add time
    %----------------------------------------------------------------------
    
    fixPlot(hRa,rMin,rMax,RMin,RMax,rText,rhoText,plotTime,legPos,legTextR);
    
    if(~isempty(legTextP))  
        fixPlot(hPa,pMin,pMax,PMin,PMax,rText,vText,plotTime,legPos,legTextP);
    end
    
    % write the figure files
    save2pdf(outputFile,hRPf);
    %close(hRPf);

end % for iPlot

fprintf(1,'Finished\n');
