function plotEquilibria2D(stoc,ddft,optsPlotGIF,equilibria,pdfFile)
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
doP=~isempty(optsPlotGIF(1).legTextP);

% set initial and final time positions
plotPos(2)=1;
plotPos(3)=length(optsPlotGIF(1).plotTimes);

fullscreen = get(0,'ScreenSize');

for iPlot=2:2  % we're using optsPlotGIF(2) and optsPlotGIF(3)

    % choose appropriate plotting options
    optsPlot=optsPlotGIF(iPlot);

    %----------------------------------------------------------------------
    % Set up figure
    %----------------------------------------------------------------------
    
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

%     if(~strcmp(optsPlot.oneLeg,'off'))
%         % create the legend
%         addOneLeg(optsPlot,hRPf);
%     end
  
    %----------------------------------------------------------------------
    % Get time, colours and types
    %----------------------------------------------------------------------
    
    % get time to add to plot
    plotTimes=optsPlot.plotTimes;
    plotTime=plotTimes(plotPos(iPlot));

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
    outputFile=pdfFile{iPlot};
    outputFile=[outputFile(1:end-5) 'Eq.pdf'];
    
    %----------------------------------------------------------------------
    % Stochastic equilibrium plots
    %----------------------------------------------------------------------
     
    if(nStoc>0)
                   
        % get x and p values from structure
        if(iPlot==2)
            x=stoc(1).xInitial;
        else
            x=stoc(1).xFinal;
        end
        p=stoc(1).pIF;
            
        x=x';
        p=p';
        
        optsPlot.type=stocType(1,:);
        
        optsPlot.faceColour = optsPlot.lineColourStoc{1};
        
        if(~isempty(x))
            plotRhoVdistStoc2D(x,p,optsPlot,handlesRP(iPlot));
        end
        
        hold(hRa,'on');
        hold(hPa,'on');
        
    end
       
        
    %----------------------------------------------------------------------
    % DDFT data plots
    %----------------------------------------------------------------------
      
    if(nDDFT>0)
        
        % NEED TO OUTPUT FINAL EQUILIBRIUM
        
        if(iPlot==2)
            
            for iDDFT=1:nDDFT
            
                % get rho, v and r values
                rho=ddft(iDDFT).rho_t;

                optsPlot.type=DDFTType(1,:);

                % get values at appropriate time
                rhot=rho(:,:,1);

                if(strcmp(optsPlot.type,'rv'))
                    v=ddft(1).v_IP;
                    vt=v(:,:,plotPos(iPlot));
                else
                    v=ddft(1).flux_t;
                    fluxNorm = 0.1*max(max(max(v)));
                    optsPlot.fluxNorm=fluxNorm;

                    vt=v(:,:,plotPos(iPlot));
                end
                
                optsPlot.faceColour=optsPlot.lineColourDDFT{iDDFT};
                
                %ddft(iDDFT).shape
                
                % plot the distributions
                plotRhoVdistDDFT2D(rhot,vt,ddft(iDDFT).shape.Interp,ddft(iDDFT).shape.Pts,optsPlot,handlesRP(iPlot));

                hold(hRa,'on');
                hold(hPa,'on');
        
            end
                
        end
        
    end
            
    %----------------------------------------------------------------------
    % Set axes, legend, add time
    %----------------------------------------------------------------------
    
    optsPlot.time=plotTime;
    
    optsPlot.xLab='x';
    optsPlot.yLab='y';
    
    optsPlot.xMin=optsPlot.rMin{1};
    optsPlot.xMax=optsPlot.rMax{1};
    optsPlot.yMin=optsPlot.rMin{2};
    optsPlot.yMax=optsPlot.rMax{2};

    optsPlot.zMin=optsPlot.RMin;
    optsPlot.zMax=optsPlot.RMax;

    optsPlot.zLab='Density';

    %optsPlot.legText=optsPlot.legTextR{iSpecies};
    if(strcmp(optsPlotGIF(1).plotType,'surf'))
        fixPlot2Dsurf(hRa,optsPlot);
        optsPlot.time=[];
        fixPlot2Dcontour(hPa,optsPlot);
    else
        fixPlot2Dcontour(hRa,optsPlot);
    end
        
    if(doP)  

        optsPlot.xMin=optsPlot.pMin{1};
        optsPlot.xMax=optsPlot.pMax{1};
        optsPlot.yMin=optsPlot.pMin{2};
        optsPlot.yMax=optsPlot.pMax{2};            

        %optsPlot.legText=optsPlot.legTextP{iSpecies};
        optsPlot.zLab='Momentum';
        optsPlot.zMin=optsPlot.PMin{1};
        optsPlot.zMax=optsPlot.PMax{1};

        fixPlot2D(hPa,optsPlot);
    end
        
    
    % write the figure files
    save2pdf(outputFile,hRPf,100,true);
    %close(hRPf);

end % for iPlot

fprintf(1,'Finished\n');
