function outputFiles = plotInitialFinal2D(stoc,ddft,optsPlot,equilibria)
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

% easy way to check if any of the computations include momentum
doP=~isempty(optsPlot.legTextP);

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

outputFiles = {};

% set initial and final time positions
plotPos(1)=1;
plotPos(2)=length(optsPlot.plotTimes);

fullscreen = get(0,'screensize');

for iPlot=1:2  

    %----------------------------------------------------------------------
    % Set up figure
    %----------------------------------------------------------------------
    
    if(fullscreen(3)>1500)
        fullscreen(3)=fullscreen(3)/2;
    end
    hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
    % set background colour to white
    set(hRPf,'Color','w');

    switch optsPlot.plotType
        case 'surf'  
            handles=tightsubplot(1,2,0.075,0.075,0.075);
            hRa=handles(1);
            hPa=handles(2);
        case 'contour'
            hRa=axes;
            % invisible axes off the figure
            hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
    end

    % figure and axes handles to be passed to plotting functions
    handlesRP(iPlot)=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa); %#ok

    %----------------------------------------------------------------------
    % Set up legend
    %----------------------------------------------------------------------

    if(~strcmp(optsPlot.oneLeg,'off'))
        % create the legend
        addOneLeg(optsPlot,hRPf);
    end
  
    %----------------------------------------------------------------------
    % Get time, colours and types
    %----------------------------------------------------------------------
    
    % get time to add to plot
    plotTimes=optsPlot.plotTimes;
    plotTime=plotTimes(plotPos(iPlot));

    % surface colours and types
    if(nStoc>0)
         lineColourStoc=optsPlot.lineColourStoc;
         % get type info -- used to decide whether to plot the velocity
         stocType=optsPlot.stocType;
    else
        lineColourStoc=[];
        stocType=[];
    end
   
    if(nDDFT>0)
        lineColourDDFT=optsPlot.lineColourDDFT;
        % get type info -- used to decide whether to plot the velocity
        DDFTType=optsPlot.DDFTType;
    else
        lineColourDDFT=[];
        % get type info -- used to decide whether to plot the velocity
        DDFTType=[];
    end
    
    % and file to save in
    outputFile=optsPlot.IFFiles{iPlot};
    
    %----------------------------------------------------------------------
    % Stochastic data plots
    %----------------------------------------------------------------------
     
    for iStoc=1:nStoc
    
        optsPlot.faceColour=lineColourStoc{iStoc};
                   
        rho   = stoc(iStoc).rho(:,:,:,plotPos(iPlot));
        %v     = stoc(iStoc).v(:,:,:,:plotPos(iPlot));
        flux  = stoc(iStoc).flux(:,:,:,:,plotPos(iPlot));
        boxes =  stoc(iStoc).boxes(:,:,:,:,plotPos(iPlot));
        
        optsPlot.plotTime=plotTime;        
        optsPlot.type=stocType{iStoc};
        
        plotRhoVdistStoc2D(rho,flux,boxes,optsPlot,handlesRP(iPlot));
        
%         hold(hRa,'on');
%         hold(hPa,'on');
        
    end
       
        
    %----------------------------------------------------------------------
    % DDFT data plots
    %----------------------------------------------------------------------
      
    for iDDFT=1:nDDFT
        
        % get rho, v and r values
        rho=ddft(iDDFT).rho_t;
        
        optsPlot.type=DDFTType{iDDFT};
          
        % get values at appropriate time
        rhot=rho(:,:,plotPos(iPlot));
        
        if(strcmp(optsPlot.type,'rv'))
            v=ddft(iDDFT).v_IP;
            vt=v(:,:,plotPos(iPlot));
        else
            v=ddft(iDDFT).flux_t;
            fluxNorm = 0.1*max(max(max(v)));
            optsPlot.fluxNorm=fluxNorm;
            
            vt=v(:,:,plotPos(iPlot));
        end
        
        optsPlot.faceColour=lineColourDDFT{iDDFT};
        
        % plot the distributions
        plotRhoVdistDDFT2D(rhot,vt,ddft(iDDFT).shape.Interp,ddft(iDDFT).shape.Pts,optsPlot,handlesRP(iPlot));

        hold(hRa,'on');
        hold(hPa,'on');
    end
            
    %----------------------------------------------------------------------
    % Set axes, legend, add time
    %----------------------------------------------------------------------
    
    optsPlot.time=plotTime;
    
    optsPlot.xLab='x';
    optsPlot.yLab='y';
    
    optsPlot.xMin=optsPlot.rMin(1);
    optsPlot.xMax=optsPlot.rMax(1);
    optsPlot.yMin=optsPlot.rMin(2);
    optsPlot.yMax=optsPlot.rMax(2);

    optsPlot.zMin=optsPlot.RMin;
    optsPlot.zMax=optsPlot.RMax;

    optsPlot.zLab='Density';

    %optsPlot.legText=optsPlot.legTextR{iSpecies};
    if(strcmp(optsPlot.plotType,'surf'))
        fixPlot2Dsurf(hRa,optsPlot);
        optsPlot.time=[];
        fixPlot2Dcontour(hPa,optsPlot);
    else
        fixPlot2Dcontour(hRa,optsPlot);
    end
        
    if(doP)  

        optsPlot.xMin=optsPlot.pMin(1);
        optsPlot.xMax=optsPlot.pMax(1);
        optsPlot.yMin=optsPlot.pMin(2);
        optsPlot.yMax=optsPlot.pMax(2);            

        %optsPlot.legText=optsPlot.legTextP{iSpecies};
        optsPlot.zLab='Momentum';
        optsPlot.zMin=optsPlot.PMin(1);
        optsPlot.zMax=optsPlot.PMax(1);

        fixPlot2D(hPa,optsPlot);
    end
        
    
    % write the figure files
    save2pdf(outputFile,hRPf,100,true);
    %close(hRPf);
    
    outputFiles = cat(2,outputFiles,outputFile);

end % for iPlot

fprintf(1,'Finished\n');
