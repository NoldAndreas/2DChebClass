function outputFile = plotEquilibria2D(stoc,ddft,optsPlot,optsPhys,equilibria)
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
doP=~isempty(optsPlot(1).legTextP);

fullscreen = get(0,'ScreenSize');
%----------------------------------------------------------------------
% Set up figure
%----------------------------------------------------------------------

hRPf=figure('Position',[0 0 0.75*fullscreen(3) 0.75*fullscreen(4)]);
% set background colour to white
set(hRPf,'Color','w');


separateSpecies = optsPlot.separateSpecies;

separateError = optsPlot.separateError;


if(separateSpecies)
    nSpecies = length(optsPhys.nParticlesS);
end

switch optsPlot.plotType
    case 'surf'  
        if(~separateSpecies)
            handles=tightsubplot(1,2,0.075,0.075,0.075);
            hRa=handles(1);
            hPa=handles(2);
            % figure and axes handles to be passed to plotting functions
            handlesRP = struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa);
        else
            if(~separateError)
                handles=tightsubplot(2,nSpecies,0.075,0.075,0.075);
                hRa = handles(1:nSpecies);
                hPa = handles(nSpecies+1:end);
                for iSpecies = 1:nSpecies
                    handlesRP(iSpecies) = struct('hRPf',hRPf,'hRa',hRa(iSpecies),'hPa',hPa(iSpecies));%#ok
                end
            else
                hETemp = axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
                
                handles=tightsubplot(3,nSpecies,0.075,0.075,0.075);
                hRa = handles(1:nSpecies);
                hEa = handles(nSpecies+1:2*nSpecies);
                hPa = handles(2*nSpecies+1:end);
                for iSpecies = 1:nSpecies
                    handlesRP(iSpecies) = struct('hRPf',hRPf,'hRa',hRa(iSpecies),'hPa',hPa(iSpecies));%#ok
                    handlesRPE(iSpecies) = struct('hRPf',hRPf,'hRa',hEa(iSpecies),'hPa',hETemp);%#ok
                end
            end

        end
    case 'contour'
        if(~separateSpecies)
            hRa=axes;
            % invisible axes off the figure
            hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
        else
            hRa = tightsubplot(1,nSpecies,0.075,0.075,0.075);
            %hPa = hPa*ones(1,nSpecies);
            for iSpecies = 1:nSpecies
                hPa(iSpecies) = axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');%#ok
                handlesRP(iSpecies) = struct('hRPf',hRPf,'hRa',hRa(iSpecies),'hPa',hPa(iSpecies));%#ok
            end
        end          
end

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
plotTime=plotTimes(1);

% surface colours and types
if(nStoc>0)
     optsPlot.faceColour=optsPlot.eqColour;
     % get type info -- used to decide whether to plot the velocity
     stocType=optsPlot.stocType;
end

if(nDDFT>0)
    optsPlot.faceColour=optsPlot.eqColour;
    % get type info -- used to decide whether to plot the velocity
    DDFTType=optsPlot.DDFTType;
end

% and file to save in
outputFile = optsPlot.eqFile;


%----------------------------------------------------------------------
% DDFT data plots
%----------------------------------------------------------------------

if(nDDFT>0)

    for iDDFT = 1:nDDFT
    
        % get rho, v and r values
        rho = ddft(iDDFT).dynamicsResult.rho_t;
        flux = ddft(iDDFT).dynamicsResult.flux_t;

        optsPlot.type=DDFTType(iDDFT);
        optsPlot.lineStyle = optsPlot.lineStyleDDFT{iDDFT};

        % get values at appropriate time
        rhot  = rho(:,:,1);
        fluxt = flux(:,:,:,1);

        optsPlot.fluxNorm = 1;

        colours = optsPlot.lineColourDDFT{iDDFT};

        if(~separateSpecies)
            optsPlot.faceColour = colours;
            % plot the distributions
            plotRhoVdistDDFT2D(rhot,fluxt,ddft(iDDFT).IDC,optsPlot,handlesRP);

            hold(hRa,'on');
            hold(hPa,'on');
        else
            for iSpecies = 1:nSpecies
                optsPlot.faceColour = colours(iSpecies);
                optsPlot.doFlux = false;
                plotRhoVdistDDFT2D(rhot(:,iSpecies),fluxt(:,iSpecies),ddft(iDDFT).IDC,optsPlot,handlesRP(iSpecies));
                hold(hRa(iSpecies),'on');
                hold(hPa(iSpecies),'on');
            end
        end
    end

end

%----------------------------------------------------------------------
% Stochastic equilibrium plots
%----------------------------------------------------------------------

if(nStoc>0)

    %for iStoc = 1:nStoc
    for iStoc = 1:1   % only one stochastic equilibrium
    
        rho   = equilibria(iStoc).data.REq;
        v     = equilibria(iStoc).data.vEq;
        boxes = equilibria(iStoc).data.xEq;

        if(separateError)
            y1C = boxes(:,:,1,1);
            y2C = boxes(:,:,2,1);
            y1C_kv = y1C(:);
            y2C_kv = y2C(:);
            IDC = ddft(iDDFT).IDC;

            pts         = IDC.GetInvCartPts(y1C_kv,y2C_kv);

            InterPol = IDC.InterpolationMatrix_Pointwise(pts.y1_kv,pts.y2_kv);

            % only compare first DDFT - to generalise
            rhoDDFT = ddft(1).dynamicsResult.rho_t;
            rhoDDFT = rhoDDFT(:,:,1);
            
            rhoDDFTBoxes = InterPol*rhoDDFT;
        
        
            nBins = optsPlot.nBins;
            N1 = nBins(1);
            N2 = nBins(2);

%             rho3 = rhoIBoxes(:,3);
%             rho3 = reshape(rho3,N1,N2);
%         
%         
%             surf(hEa(3),y1C,y2C,rho3)
        
        end
       
        
        optsPlot.type=stocType(iStoc);
        optsPlot.lineStyle = optsPlot.lineStyleStoc{iStoc};

        colours = optsPlot.lineColourStoc{iStoc};
        nParticlesS = optsPlot.nParticlesS;

        if(~separateSpecies)
            optsPlot.faceColour = colours;
            plotRhoVdistStoc2D(rho,v,boxes,optsPlot,handlesRP);
            hold(hRa,'on');
            hold(hPa,'on');
        else
            if(~separateError)
            
                for iSpecies = 1:nSpecies
                    optsPlot.faceColour = colours(iSpecies);
                    optsPlot.nParticlesS = nParticlesS(iSpecies);
                    plotRhoVdistStoc2D(rho(:,:,iSpecies),v,boxes,optsPlot,handlesRP(iSpecies));
                    hold(hRa(iSpecies),'on');
                    hold(hPa(iSpecies),'on');
                end
                
            else
                
                maxError = zeros(nSpecies,1);
                
                for iSpecies = 1:nSpecies
                    rhoDDFTBoxesS = reshape( rhoDDFTBoxes(:,iSpecies), N1,N2);
                    rhoError = abs(rho(:,:,iSpecies) - rhoDDFTBoxesS);
                    maxError(iSpecies) = max(max(rhoError));
                    optsPlot.faceColour = colours(iSpecies);
                    optsPlot.nParticlesS = nParticlesS(iSpecies);
                    plotRhoVdistStoc2D(rhoError,v,boxes,optsPlot,handlesRPE(iSpecies));
                    hold(hEa(iSpecies),'on');
                    hold(hPa(iSpecies),'on');
                    
                end
            end
                
        end
        
    end

end


%----------------------------------------------------------------------
% Set axes, legend, add time
%----------------------------------------------------------------------

%optsPlot.time=plotTime;
optsPlot.time=[];

optsPlot.xLab='$y_1$';
optsPlot.yLab='$y_2$';

optsPlot.xMin=optsPlot.rMin(1);
optsPlot.xMax=optsPlot.rMax(1);
optsPlot.yMin=optsPlot.rMin(2);
optsPlot.yMax=optsPlot.rMax(2);

optsPlot.zMin=optsPlot.RMin;
optsPlot.zMax=optsPlot.RMax;

optsPlot.zLab='Density';

%optsPlot.legText=optsPlot.legTextR{iSpecies};
if(strcmp(optsPlot(1).plotType,'surf'))
    if(~separateSpecies)
        fixPlot2Dsurf(hRa,optsPlot);
        optsPlot.time=[];
        fixPlot2Dcontour(hPa,optsPlot);
    else
        for iSpecies = 1:nSpecies
            fixPlot2Dsurf(hRa(iSpecies),optsPlot);
            
            if(separateError)
                optsPlotError = optsPlot;
                optsPlotError.zLab = 'Error';
                optsPlotError.zMax = 1.1*max(maxError);
                fixPlot2Dsurf(hEa(iSpecies),optsPlotError);
                
            end
            
            optsPlot.time=[];
            fixPlot2Dcontour(hPa(iSpecies),optsPlot);
        end
    end
else
    
    if(~separateSpecies)
        fixPlot2Dcontour(hRa,optsPlot);
    else
        for iSpecies = 1:nSpecies
            fixPlot2Dcontour(hRa(iSpecies),optsPlot);
        end
    end

    

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

    if(~separateSpecies)
        fixPlot2D(hPa,optsPlot);
    else
        for iSpecies = 1:nSpecies
            fixPlot2D(hPa(iSpecies),optsPlot);
        end
    end
    
end


% write the figure files
save2pdf(outputFile,hRPf,100,true);
%close(hRPf);


fprintf(1,'Finished\n');
