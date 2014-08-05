function outputFile = makeMoviePlanar2D(stoc,ddft,optsPlot,optsPhys,equilibria)
% makeMovie(stoc,ddft,optsPlot,equilibria)
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

% get times
plotTimes=optsPlot.plotTimes;

if(nStoc>0)
    % choose colours
    lineColourStoc=optsPlot.lineColourStoc;
    % get type info -- used to decide whether to plot the velocity
    stocType=optsPlot.stocType;
else
    lineColourStoc=[];
    stocType=[];
end


if(nDDFT>0)
    % choose colours
    lineColourDDFT=optsPlot.lineColourDDFT;
    % get type info -- used to decide whether to plot the velocity
    DDFTType=optsPlot.DDFTType;
else
 lineColourDDFT=[];

 DDFTType=[];
end


separateSpecies = optsPlot.separateSpecies;
separateComp = optsPlot.separateComp;

% determine whether to do any momentum plotting
doP=~isempty(optsPlot.legTextP);

% number of plots in the movie
nPlots=length(plotTimes);

% determine which movies to make
doGif=optsPlot.doMovieGif;
doSwf=optsPlot.doMovieSwf;
doPdfs=optsPlot.doPdfs;

% get movie files and options
pdfDir=optsPlot.pdfDir;
movieFile=[optsPlot.movieFile optsPlot.plotType ...
           '-' num2str(separateSpecies) '-' num2str(separateComp)];
fps=optsPlot.fps;
dpi=optsPlot.dpi;


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

%fullscreen = get(0,'ScreenSize');

%hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
hRPf=figure('Position',[0 0 1200 800]);
% set background colour to white
set(hRPf,'Color','w');




nSpecies = length(optsPhys.nParticlesS);
nComp = max(nStoc,nDDFT);

if(separateSpecies)
    nAxesS = nSpecies;
else
    nAxesS = 1;
end
if(separateComp)
    nAxesC = nComp;
else
    nAxesC = 1;
end

nAxes = nAxesS*nAxesC;

switch optsPlot.plotType
    case 'surf'  
            handles=tightsubplot(nAxes,2,0.075,0.075,0.075);
            hRa = handles(1:2:end);
            hPa = handles(2:2:end);
        
    case 'contour'
            hPa = axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
            hPa = hPa*ones(1,nAxes);
end

if(separateSpecies && separateComp)
    hRaTemp = reshape(hRa,nAxesC,nAxesS);
    hPaTemp = reshape(hPa,nAxesC,nAxesS);
elseif(separateSpecies)
    hRaTemp = zeros(nComp,nSpecies);
    hPaTemp = zeros(nComp,nSpecies);
    for iComp = 1:nComp
        hRaTemp(iComp,:) = hRa;
        hPaTemp(iComp,:) = hPa;
    end
elseif(separateComp)
    hRaTemp = zeros(nComp,nSpecies);
    hPaTemp = zeros(nComp,nSpecies);
    for iSpecies = 1:nSpecies
        hRaTemp(:,iSpecies) = hRa.';
        hPaTemp(:,iSpecies) = hPa.';
    end    
else
    hRaTemp = hRa*ones(nComp,nSpecies);
    hPaTemp = hPa*ones(nComp,nSpecies);
end
for iComp = 1:nComp
    handlesRP(iComp)=struct('hRPf',hRPf,'hRa',hRaTemp(iComp,:),'hPa',hPaTemp(iComp,:)); %#ok
end
        
set(hRa,'nextplot','replacechildren');
set(hPa,'nextplot','replacechildren');


if(doGif)
    gifFile=[movieFile '.gif'];
    delayTime=1/fps;
end

%--------------------------------------------------------------------------
% Set up legend
%--------------------------------------------------------------------------

if(~strcmp(optsPlot.oneLeg,'off'))
    % create the legend
    addOneLeg(optsPlot,hRPf);
end


%--------------------------------------------------------------------------
% Make movie frame
%--------------------------------------------------------------------------


for iPlot=1:nPlots

    for iAxis = 1:nAxes
        hold(hRa(iAxis),'off');
        hold(hPa(iAxis),'off')
    end

    % get plot time
    plotTime=plotTimes(iPlot);
    optsPlot.plotTime=plotTime;

    %----------------------------------------------------------------------
    % Stochastic plotting
    %----------------------------------------------------------------------
    
    for iStoc=1:nStoc
           
        optsPlot.faceColour=lineColourStoc{iStoc};   
       
        rho   = stoc(iStoc).rho(:,:,:,iPlot);
        flux  = stoc(iStoc).flux(:,:,:,:,iPlot);
        boxes =  stoc(iStoc).boxes(:,:,:,:,iPlot);
        
        optsPlot.plotTime=plotTime;        
        optsPlot.type=stocType{iStoc};
        
        plotRhoVdistStoc2D(rho,flux,boxes,optsPlot,handlesRP(iStoc));
            
%         hold(hMRa,'on');
%         hold(hMPa,'on');
        
    end

    %----------------------------------------------------------------------
    % DDFT plotting
    %----------------------------------------------------------------------
    
    for iDDFT=1:nDDFT

        % get rho, v and r values
        rho=ddft(iDDFT).dynamicsResult.rho_t;

        optsPlot.type=DDFTType{iDDFT};
        
        optsPlot.faceColour=lineColourDDFT{iDDFT};
           
        % get values at appropriate time
        rhot=rho(:,:,iPlot);
        
        if(strcmp(optsPlot.type,'rv'))
            v=ddft(iDDFT).dynamicsResult.v_IP;
            vt=v(:,:,iPlot);
        else
            v=ddft(iDDFT).dynamicsResult.flux_t;
            %fluxNorm = 0.1*max(max(max(abs(v))));
            v1=v(1:end/2,:,:);
            v2=v(end/2+1:end,:,:);
            v1=v1(:);
            v2=v2(:);
            vNorm=sqrt(v1.^2+v2.^2);
            fluxNorm = max(vNorm);
            optsPlot.fluxNorm=fluxNorm;
            
            vt=v(:,:,iPlot);
        end
            
        % plot the distributions
        plotRhoVdistDDFT2D(rhot,vt,ddft(iDDFT).IDC.Interp,ddft(iDDFT).IDC.Pts,optsPlot,handlesRP(iDDFT));
        
%         for iAxis = 1:nAxes
%             hold(hRa(iAxis),'on');
%             hold(hPa(iAxis),'on')
%         end
    
    end
        
    %----------------------------------------------------------------------
    % Fix axes, labels, etc
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

    % FIX THIS FOR STOCHASTIC
 
    for iAxis = 1:nAxes
    
        switch optsPlot.plotType
            case 'surf'
                fixPlot2Dsurf(hRa(iAxis),optsPlot);
                optsPlot.time=[];
                fixPlot2Dcontour(hPa(iAxis),optsPlot);

                optsPlot.time=plotTime;
                
                if(nDDFT>0)
                    optsPlot.lineColour=lineColourDDFT{1};
                else
                    optsPlot.lineColour=lineColourStoc{1};
                end
                if(iPlot==1)
                   % optsPlot.quiverNorm=getQuiverNorm(optsPlot);
                   otpsPlot.quiverNorm = 1;
                end

%                addPotential(hRa(iAxis),optsPlot);

            case 'contour'
                fixPlot2Dcontour(hRa(iAxis),optsPlot);
                
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
        
        fixPlot2D(hPa,optsPlot);

    end
 
    %----------------------------------------------------------------------
    % Get gif frames and save write file
    %----------------------------------------------------------------------
    
    if(doGif)
        
        % plot the current frame
        f = getframe(hRPf);

        % get image data
        % need new map each frame as colours can change
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        
        if(iPlot==1)
            %make poster file for beamer presentations
            posterFile=[movieFile 'Poster.png'];
            imwrite(im,map,posterFile,'png');
            
            % write first frame
            imwrite(im,map,gifFile,'DelayTime',delayTime,'LoopCount',inf);
        else
            % append this frame
            [im,map] = rgb2ind(f.cdata,256,'nodither');  
            imwrite(im,map,gifFile,'gif','WriteMode','append','DelayTime',delayTime);
        end
        
        outputFile = gifFile;
        
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
    
    outputFile = swfFile;
    
end

close(hRPf);