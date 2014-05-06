function plotInitialFinal2D(stoc,ddft,optsPlotGIF,equilibria,pdfFile)
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
doP=~isempty(optsPlotGIF(1).legTextP);

% set initial and final time positions
plotPos(2)=1;
plotPos(3)=length(optsPlotGIF(1).plotTimes);

handlesRP=struct('hRPf',0,'hRa',0,'hPa1',0,'hPa2',0);

for iPlot=2:3  % we're using optsPlotGIF(2) and optsPlotGIF(3)

    % choose appropriate plotting options
    optsPlot=optsPlotGIF(iPlot);

    % get time to add to plot
    plotTimes=optsPlot.plotTimes;
    plotTime=plotTimes(plotPos(iPlot));

%     % line styles and types
    if(nStoc>0)
%         lineStyleStoc=optsPlot.lineStyleStoc;
%         lineMarkerStoc=optsPlot.lineMarkerStoc;
         lineColourStoc=optsPlot.lineColourStoc;
%         % get type info -- used to decide whether to plot the velocity
         stocType=optsPlot.stocType;
    else
%         % choose line styles
%         lineStyleStoc=[];
%         lineMarkerStoc=[];
%         lineColourStoc=[];
%         % get type info -- used to decide whether to plot the velocity
         stocType=[];
     end
%     
     if(nDDFT>0)
%         lineStyleDDFT=optsPlot.lineStyleDDFT;
%         lineMarkerDDFT=optsPlot.lineMarkerDDFT;
         lineColourDDFT=optsPlot.lineColourDDFT;
%          lineColourDDFT={{'b'}};
%         % get type info -- used to decide whether to plot the velocity
         DDFTType=optsPlot.DDFTType;
%         DDFTType='r';
         
     else
%         lineStyleDDFT=[];
%         lineMarkerDDFT=[];
         lineColourDDFT=[];
%         % get type info -- used to decide whether to plot the velocity
         DDFTType=[];
     end

%     % get legend settings
%     legPos=optsPlot.legPos;
%     legTextR=optsPlot.legTextR;
%     legTextP=optsPlot.legTextP;
    
    % and file to save in
    outputFile=pdfFile{iPlot};
    
    % set up figure
    fullscreen = get(0,'ScreenSize');
    hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
    % set background colour to white
    set(hRPf,'Color','w');
    
    % if any of the calculations contain momentum, plot them to a visible axis,
    % otherwise create a dummy set of axes which aren't visible
    
    nSpecies=length(optsPlot.nParticlesS);
    
    hRa=zeros(nSpecies);
    hPa1=zeros(nSpecies);
    hPa2=zeros(nSpecies);
    
    if(doP)    
        for iSpecies=1:nSpecies
            hRa(iSpecies)=subplot(nSpecies,3,1+3*(iSpecies-1));
            hPa1(iSpecies)=subplot(nSpecies,3,2+3*(iSpecies-1));
            hPa2(iSpecies)=subplot(nSpecies,3,3+3*(iSpecies-1));
        end
    else     
        % invisible axes off the figure
        hPaTemp= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
       
        for iSpecies=1:nSpecies
            hRa(iSpecies)=subplot(1,nSpecies,iSpecies);
            hPa1(iSpecies)=hPaTemp;
            hPa2(iSpecies)=hPaTemp;
        end

    end
     
    % figure and axes handles to be passed to plotting functions
    handlesRP(iPlot)=struct('hRPf',hRPf,'hRa',hRa,'hPa1',hPa1,'hPa2',hPa2);

    %----------------------------------------------------------------------
    % Legend
    %----------------------------------------------------------------------
    
    % Think about this, but we probably don't want one
    % as it's easier to do it with titles
    
%     % if we're adding one legend, do it now
%     if(~strcmp(optsPlot.oneLeg,'off'))
%         % combine legend labels, styles and colours, with the assumption
%         % that they are all contained in legTextR (i.e. all are plotted in
%         % the position plot
%         
%         optsPlot.legText=optsPlot.legTextR;
%         
%         optsPlot.lineStyle=cat(2,lineStyleStoc,lineStyleDDFT);
%         optsPlot.lineMarker=cat(2,lineMarkerStoc,lineMarkerDDFT);
%         optsPlot.lineColour=cat(2,lineColourStoc,lineColourDDFT);
%         
%         % create the legend
%         addOneLeg(optsPlot,hRPf);
%     end 
    
    %----------------------------------------------------------------------
    % Stochastic data plots
    %----------------------------------------------------------------------
     
    for iStoc=1:nStoc
    
        optsPlot.faceColour=lineColourStoc{iStoc};
               
%         optsPlot.lineStyle=lineStyleStoc{iStoc};
%         optsPlot.lineMarker=lineMarkerStoc{iStoc};
        
        % get x and p values from structure
        x=stoc(iStoc).x;
        p=stoc(iStoc).p;
        
        % get values at appropriate time
        % transpose as should be nParticles x nSamples
        xt=x(:,:,plotPos(iPlot))'; 
        pt=p(:,:,plotPos(iPlot))';
        
        optsPlot.type=stocType(iStoc,:);
        
        plotRhoVdistStoc2D(xt,pt,optsPlot,handlesRP(iPlot));
        
        hold(hRa,'on');
        hold(hPa1,'on');
        hold(hPa2,'on');
        
    end
       
        
    %----------------------------------------------------------------------
    % DDFT data plots
    %----------------------------------------------------------------------
    
    % FIX THIS
    
    for iDDFT=1:nDDFT
        
        % get rho, v and r values
        rho=ddft(iDDFT).rho_t;
        
        optsPlot.type=DDFTType(iStoc,:);
        
        % NEED TO CHANGE THIS WHEN WE HAVE MULTIPLE SPECIES
        
        % get values at appropriate time
        %rhot=rho(:,:,plotPos(iPlot));
        rhot=rho(:,plotPos(iPlot));
        
        if(strcmp(optsPlot.type,'rv'))
            v=ddft(iDDFT).v_IP;
            %vt=v(:,:,plotPos(iPlot));
            vt=v(:,plotPos(iPlot));
        else
            vt=[];
        end
        
        %         % choose line style
%         optsPlot.lineStyle=lineStyleDDFT{iDDFT};
%         optsPlot.lineMarker=lineMarkerDDFT{iDDFT};
        optsPlot.faceColour=lineColourDDFT{iDDFT};
        
        % plot the distributions
        plotRhoVdistDDFT2D(rhot,vt,ddft(iDDFT).Interp,optsPlot,handlesRP(iPlot));

        hold(hRa(iSpecies),'on');
        hold(hPa1(iSpecies),'on');
        hold(hPa2(iSpecies),'on');

    end
    
%     if(doEqDDFT(iPlot))
% 
%         if(iPlot==2)
%             doI=true;
%             doF=false;
%         else
%             doI=false;
%             doF=true;
%         end
%         
%         addEqLinesDDFT(handlesRP,[],ddft,optsPlotGIF(iPlot),doI,doF)
% 
%     end
    
    %----------------------------------------------------------------------
    % Equilibrium stochastic data plots
    %----------------------------------------------------------------------
    
    if(nStoc>0 && optsPlotGIF(1).plotEquilibria)
        
        optsPlot.type=stocType(1,:);
        optsPlot.faceColour=optsPlot.eqColour{1};
             
        pE=stoc(1).pIF;
        
        if(iPlot==2)
            xI=stoc(1).xInitial;
            plotRhoVdistStoc2D(xI',pE',optsPlot,handlesRP(2)); 
        end
        
        if(iPlot==3)
            xF=stoc(1).xFinal;
            if(~isempty(xF))
                plotRhoVdistStoc2D(xF',pE',optsPlot,handlesRP(3));
            end
        end
     
    end


    %----------------------------------------------------------------------
    % Set axes, legend, add time
    %----------------------------------------------------------------------
    
    %fixPlot(hRa,rMin,rMax,RMin,RMax,rText,rhoText,plotTime,legPos,legTextR);
    
    optsPlot.time=plotTime;
    
    % FIX THIS
    %optsPlot.viewPoint=[31.5 18];
    optsPlot.xLab='x';
    optsPlot.yLab='y';
    
    optsPlot.xMin=optsPlot.rMin{1};
    optsPlot.xMax=optsPlot.rMax{1};
    optsPlot.yMin=optsPlot.rMin{2};
    optsPlot.yMax=optsPlot.rMax{2};

    optsPlot.zMin=optsPlot.RMin;
    optsPlot.zMax=optsPlot.RMax;


    for iSpecies=1:nSpecies

        optsPlot.zLab=['Density; Species ' num2str(iSpecies) '  '];
        
        %optsPlot.legText=optsPlot.legTextR{iSpecies};
        if(strcmp(optsPlotGIF(1).plotType,'surf'))
            fixPlot2D(hRa(iSpecies),optsPlot);
        else
            fixPlot2Dcontour(hRa(iSpecies),optsPlot);
        end
    end

    if(doP)  

        optsPlot.xMin=optsPlot.pMin{1};
        optsPlot.xMax=optsPlot.pMax{1};
        optsPlot.yMin=optsPlot.pMin{2};
        optsPlot.yMax=optsPlot.pMax{2};            

        %optsPlot.legText=optsPlot.legTextP{iSpecies};
        optsPlot.zLab=['Momentum 1; Species ' num2str(iSpecies) '  '];
        optsPlot.zMin=optsPlot.PMin{1};
        optsPlot.zMax=optsPlot.PMax{1}; 

        for iSpecies=1:nSpecies
            fixPlot2D(hPa1(iSpecies),optsPlot);
        end

        optsPlot.zLab=['Momentum 2; Species ' num2str(iSpecies) '  '];
        optsPlot.zMin=optsPlot.PMin{2};
        optsPlot.zMax=optsPlot.PMax{2}; 

        for iSpecies=1:nSpecies
            fixPlot2D(hPa2(iSpecies),optsPlot);
        end


    end
        
    
    % write the figure files
    %save2pdf(outputFile,hRPf,100,true);
    save2pdf(outputFile,hRPf);
    %close(hRPf);

end % for iPlot

fprintf(1,'Finished\n');
