function makeMovie1D(stoc,ddft,optsPlot,equilibria)
% makeMovie(stoc,ddft,optsPlotGIF,equilibria)
%   makes movie from given stochastic and DDFT data
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
%             doMovieGif      (true/false whether to make gif)
%             doPdfs          (true/false whether to save pdfs for swf)
%             doMovieSwf      (true/false whether to make swf)
%             pdfDir          (string for directory in which to save pdfs)
%             movieFile       (string for output file name (NO EXTENSION) )
%             fps             (integer frames per second)
%             dpi             (integer resolution)
%             bitmap          (true/false whether to flatten pdfs -- needed
%                              for large animations)
%             quiet           (true/false whether to suppress output of gs
%                              and pdf2swf)
%
%  equilbria -- structure of size [1 1] containing
%             REq             (vector of rho at positions xR)
%             vEq             (vector of v at positions xR)
%             xR              (vector positions, centres of histogram bins)
%             rMeanEq         (vector of mean positions, 
%                               same length as plotTimes)
%             vMeanEq         (vector of mean momenta, 
%                               same length as plotTimes)
%             plotTimes       (vector of plot times)

fprintf(1,'Making movie plots ... ');

%--------------------------------------------------------------------------
% Get plotting options
%--------------------------------------------------------------------------

% number of stochastic and DDFT plots to do
nStoc=size(stoc,2);
nDDFT=size(ddft,2);

% species masses
mS=optsPlot.mS;

% get times
plotTimes=optsPlot.plotTimes;

% get initial and final times
tMin=plotTimes(1);
tMax=plotTimes(end);


if(nStoc>0)
    % choose line styles
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
    % choose line styles
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

% determine whether to do any momentum plotting
doP=~isempty(legTextP);

% number of plots in the movie
nPlots=length(plotTimes);

% determine whether we are plotting final equilibrium values
if(~isempty(equilibria))
    doEq=true;
else
    doEq=false;
end

% determine which movies to make
doGif=optsPlot.doMovieGif;
doSwf=optsPlot.doMovieSwf;
doPdfs=optsPlot.doPdfs;

% get movie files and options
plotDir = optsPlot.plotDir;
pdfDir=optsPlot.pdfDir;
movieFile=optsPlot.movieFile;
fps=optsPlot.fps;
dpi=optsPlot.dpi;

if(~exist(optsPlot.plotDir,'dir'))
    mkdir(plotDir);
end   

% set up pdf file names
if(doPdfs || doSwf)
    if(~exist(pdfDir,'dir'))
        mkdir(pdfDir);
    end       
    pdfFileNames=[];
    nDigits=ceil(log10(nPlots));
    nd=['%0' num2str(nDigits) 'd'];
end

%--------------------------------------------------------------------------
% Set up figure
%--------------------------------------------------------------------------

% set up figure
fullscreen = get(0,'ScreenSize');
hRPf=figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
% set background colour to white
set(hRPf,'Color','w');

if(doP)    
    % tightsubplot wastes less space than subplot
    %handles=tightsubplot(2,2,0.075,0.06,0.05);
    %handles=tightsubplot(2,2,0.075,0.08,0.06);
    handles=tightsubplot(2,2,0.09,0.08,0.08);
    hRa=handles(1);
    hPa=handles(3);
    hMRa=handles(2);
    hMPa=handles(4);
else
    % tightsubplot wastes less space than subplot
    handles=tightsubplot(1,2,0.075,0.1,0.05);
    hRa=handles(1);
    hMRa=handles(2);
    hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
    hMPa= hPa;
end

% figure and axis handles to be passed to plotting functions
handlesRP=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa);
handlesM=struct('hRPf',hRPf,'hRa',hMRa','hPa',hMPa);

% stops flashing
set(hRa,'nextplot','replacechildren');
set(hPa,'nextplot','replacechildren');
set(hMRa,'nextplot','replacechildren');
set(hMPa,'nextplot','replacechildren');

% if we're adding one legend, do it now
if(~strcmp(optsPlot.oneLeg,'off'))
    % create the legend
    addOneLeg(optsPlot,hRPf);
end

%--------------------------------------------------------------------------
% Get axis limits
%--------------------------------------------------------------------------

% x axis limits
rMin=optsPlot.rMin;
rMax=optsPlot.rMax;
% y axis limits for distribution plots
RMin=optsPlot.RMin;
RMax=optsPlot.RMax;
% y axis limits for mean plots
RMMin=optsPlot.RMMin;
RMMax=optsPlot.RMMax;

if(doP) 
    % x axis limits
    pMin=optsPlot.pMin;
    pMax=optsPlot.pMax;
    % y axis limits for distribution plots
    PMin=optsPlot.PMin;
    PMax=optsPlot.PMax;
    % y axis limits for mean plots
    PMMin=optsPlot.PMMin;
    PMMax=optsPlot.PMMax;
end

%--------------------------------------------------------------------------
% Make movie frame
%--------------------------------------------------------------------------

nParticlesS = optsPlot.nParticlesS;
nParticles  = sum(nParticlesS);
nSpecies    = length(nParticlesS);
dim         = optsPlot.dim;


speciesMask=false(nParticles*dim,nSpecies);

speciesStart=[1; cumsum(nParticlesS)*dim+1];
speciesEnd=speciesStart(2:end)-1;

for iSpecies=1:nSpecies
    speciesMask(speciesStart(iSpecies):speciesEnd(iSpecies),iSpecies) = true;
end

for iPlot=1:nPlots
    
    hold(hRa,'off');
    hold(hPa,'off');
    hold(hMRa,'off');
    hold(hMPa,'off');

    % get plot time
    plotTime=plotTimes(iPlot);
   
    %----------------------------------------------------------------------
    % Stochastic plotting
    %----------------------------------------------------------------------
    
    for iStoc=1:nStoc
        
        optsPlot.lineColour=lineColourStoc{iStoc};
        optsPlot.lineStyle=lineStyleStoc{iStoc};
        optsPlot.lineMarker=lineMarkerStoc{iStoc};
        
        plotRVmeans(stoc(iStoc).rMean,stoc(iStoc).fluxMean,plotTimes,iPlot,optsPlot,handlesM,stocType(iStoc));

        hold(hMRa,'on');
        hold(hMPa,'on');        
        
        optsPlot.plotTime=plotTime;
        
        rho   = stoc(iStoc).rho(:,:,iPlot);
        flux  = stoc(iStoc).flux(:,:,iPlot);
        boxes = stoc(iStoc).boxes(:,:,iPlot);
        
        plotRhoVdistStoc(rho,flux,boxes,optsPlot,handlesRP,stocType(iStoc));
                      
    end

    %----------------------------------------------------------------------
    % DDFT plotting
    %----------------------------------------------------------------------
    
    for iDDFT=1:nDDFT
        
        % get rho, v, r and w (for means) values
        rho_t  = ddft(iDDFT).rho_t;
        flux_t = ddft(iDDFT).flux_t;
        
        Interp=ddft(iDDFT).shape.Interp;
        
        % choose line style
        optsPlot.lineStyle=lineStyleDDFT{iDDFT};
        optsPlot.lineMarker=lineMarkerDDFT{iDDFT};
        optsPlot.lineColour=lineColourDDFT{iDDFT};
        
        % plot means of r and v up to current time
        plotRVmeans(ddft(iDDFT).rMean,ddft(iDDFT).fluxMean,plotTimes,iPlot,optsPlot,handlesM,DDFTType(iDDFT));

        hold(hMRa,'on');
        hold(hMPa,'on');
        
        % get values at appropriate time
        rho  = rho_t(:,:,iPlot);
        flux = flux_t(:,:,iPlot);
        
        optsPlot.plotTime=plotTime;
        
        % plot the distributions
        %plotRhoVdistDDFT(rho,v,Interp,optsPlot,handlesRP,DDFTType(iDDFT));
        plotRhoFluxDDFT(rho,flux,Interp,optsPlot,handlesRP,DDFTType(iDDFT));
               
    end
    
    %----------------------------------------------------------------------
    % Fix axes, labels, etc
    %----------------------------------------------------------------------

    % choose appropriate text for axis labels
    if(optsPlot.symbolLabels)
        meanRText='$\bar r  \;\;$';
        meanPText='Flux';
        rhoText='$\rho \;\;$';
        vText='Flux';     
        tText='$t$';
        rText='$r$';
    else
        geom=optsPlot.geom;
        if(strcmp(geom,'spherical'))
            meanRText='Mean radial position';
            meanPText='Mean radial flux';
            vText='Radial flux';  
            rText='Radius';
        else
            meanRText='Mean position';
            meanPText='Mean flux';
            vText='Flux';  
            rText='Distance';
        end
        if(optsPlot.plotDensity)
            rhoText='Particle Density';
        else
            rhoText='Number of particles';
        end
        tText='Time';
    end
    
    fixPlot(hMRa,tMin,tMax,RMMin,RMMax,tText,meanRText,plotTime,legPos,legTextR);
    fixPlot(hRa,rMin,rMax,RMin,RMax,rText,rhoText,plotTime,legPos,legTextR);
    
    if(doP)
        fixPlot(hMPa,tMin,tMax,PMMin,PMMax,tText,meanPText,plotTime,legPos,legTextP);
        fixPlot(hPa,pMin,pMax,PMin,PMax,rText,vText,plotTime,legPos,legTextP);
    end

    %----------------------------------------------------------------------
    % Add equilibrium lines
    %----------------------------------------------------------------------
    
    if(doEq && nDDFT==0)
            for iEq = 1:size(equilibria,2)
%                addEqLines(handlesRP,handlesM,equilibria(iEq).data,optsPlot);
            end
    end
    
    if(doEq && nDDFT>0)
            ddft(1).plotTimes=plotTimes;
            addEqLinesDDFT(handlesRP,handlesM,ddft(1),optsPlot,true,true);
    end
    
    %----------------------------------------------------------------------
    % Get gif frames
    %----------------------------------------------------------------------
    
    if(doGif)
        % plot the current frame
        f = getframe(hRPf);

        % get image data
        if(iPlot==1)
            % think this makes sure that the color map is correct,
            % otherwise it is black and white
            f = getframe(hRPf);
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            
            %make poster file for beamer presentations
            posterFile=[movieFile 'Poster.png'];
            imwrite(im,map,posterFile,'png');
            
            % preallocate something in the final frame which preallocates the
            % entire matrix, giving a significant increase in performance
            im(1,1,1,nPlots) = 0;
        else
            % save image
            im(:,:,1,iPlot) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    
    %----------------------------------------------------------------------
    % Get pdf frames
    %----------------------------------------------------------------------
    
    if(doPdfs || doSwf)
        % determine output file, of the form iPlot.pdf, with leading zeros
        outputFile=[pdfDir num2str(iPlot,nd) '.pdf'];
        % add to list of pdf files used to make movie
        pdfFileNames = cat(2, pdfFileNames, [' ' outputFile]);
        
        if(doPdfs)
            % save the pdf file for this time step
            save2pdf(outputFile,hRPf,dpi);
        end
    end
        
    
end

fprintf(1,'Finished\n');

%--------------------------------------------------------------------------
% Save gif
%--------------------------------------------------------------------------

if(doGif)
    fprintf(1,'Making gif movie ... ');
    % write the movie file
    gifFile=[movieFile '.gif'];
    delayTime=1/fps;
    imwrite(im,map,gifFile,'DelayTime',delayTime,'LoopCount',inf)
    fprintf(1,'Finished\n');
end

%--------------------------------------------------------------------------
% Save swf
%--------------------------------------------------------------------------

if(doSwf)
    % name for swf file
    swfFile=[movieFile '.swf'];
    
    if(doGif)
        fprintf(1,'Making swf movie ... ');
        swfCmd=['gif2swf ' gifFile ' -o ' swfFile];
        system(swfCmd);
        fprintf(1,'Finished\n');
    else
    
        makeSwf(swfFile,pdfFileNames,optsPlot)
    end
end

%close(hRPf);