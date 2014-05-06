function handles=plotRhoVdistStoc2D(x,p,optsPlot,plotHandles)
%plotRhoVdistStoc(x,p,optsPlot,plotHandles,type)
%   plots position and velocity distributions in axes given by plotHandles
%
% INPUTS: 
%  x            -- (dim*nParticles,nSamples) matrix of positions
%  p            -- (dim*nParticles,nSamples) vector of momenta
%  optsPlot     -- a structure of size [1 1] containing:
%                        lineStyle  (linestyles for stochastic plots, 
%                                    should be  of the form ['b  ';'--b'])
%                        plotDensity (true/false whether to plot the
%                                        density; false->distruibution)
%  plotHandles  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  type         -- either 'r' for no momentum or 'rv' for momentum
%
% OUTPUTS:
%  handles      -- a structure containing the handles of the plot lines,
%                  hr for the position and if type=rv also hp for velocity

contourWidth=get(0,'defaultlinelinewidth');

nParticlesS=optsPlot.nParticlesS;
nSpecies=length(nParticlesS);

type=optsPlot.type;

if(isempty(plotHandles))
     
    fullscreen = get(0,'ScreenSize');

    hRPf=figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
    % set background colour to white
    set(hRPf,'Color','w');

    switch optsPlot.plotType
        case 'surf'
            handles=tightsubplot(2,1,0.075,0.075,0.075);
            hRa=handles(1)*ones(1,nSpecies);
            hPa=handles(2)*ones(1,nSpecies);
        case 'contour'
            hRa=axes;
            % invisible axes off the figure
            hPa= axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');
    end

    % figure and axes handles to be passed to plotting functions
    handles=struct('hRPf',hRPf,'hRa',hRa,'hPa',hPa);
      
    % otherwise assign axis handles
else
    hRa=plotHandles.hRa;
    hPa=plotHandles.hPa;
    
    if(length(hRa)<nSpecies)
        hRa = hRa(1)*ones(1,nSpecies);
    end

    if(length(hPa)<nSpecies)
        hPa = hPa(1)*ones(1,nSpecies);
    end
    
end




% determine the relevant parameters
%geom=optsPlot.geom;
nBins=optsPlot.nBins;
dim=optsPlot.dim;

% get the R and P data
%[R,P]=getRP(x,p,geom,dim);
[R,P]=getRP(x,p,'planar2D',dim);

% do histogramming of R and obtain v(r).  Note we can't just use hist as we want P=v
% as a function of r, not as a histogram

if(optsPlot.fixedBins)
    histRange=[optsPlot.rMin{1} optsPlot.rMax{1} optsPlot.rMin{2} optsPlot.rMax{2}];
    [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS,histRange);
else
    [nR,meanP,xR]= RPhist2D(R,P,nBins,nParticlesS);
end

% convert to velocities
mS=optsPlot.mS;
meanV=meanP;
for iSpecies=1:length(mS)
    meanV(:,:,:,iSpecies)=meanP(:,:,:,iSpecies)./mS(iSpecies);
end

% plotting options
faceColour=optsPlot.faceColour;

% handles of lines
hr=zeros(nSpecies,1);
% hp1=hr;
% hp2=hr;

for iSpecies=1:nSpecies
    
    nRS=nR(:,:,iSpecies);
%    meanVS=meanV(:,:,:,iSpecies);
    xRS=xR(:,:,:,iSpecies);
    
    switch optsPlot.plotType
        case 'surf'
            % surface plot
            hr(iSpecies)=surf(hRa(iSpecies),xRS(:,:,1),xRS(:,:,2),nRS);
            set(hr(iSpecies),'FaceColor',faceColour{iSpecies});
%           set(hr(iSpecies),'EdgeColor',faceColour{iSpecies});
            alpha(hr(iSpecies),0.5);
            
            hCa=hPa(iSpecies);

        case 'contour'
            
            hCa=hRa(iSpecies);
    end


    hold(hRa,'on')
    
    [C,h]=contour(hCa,xRS(:,:,1),xRS(:,:,2),nRS);
    set(h,'color',faceColour{iSpecies},'linewidth',contourWidth);
    clabel(C,h,'Color',faceColour{iSpecies});
    hold(hCa,'on');

    % FLUXES??
    
    
%     %scatter3(hRa(iSpecies),temp1,temp2,temp3,'k','filled');
%     axes(hRa(iSpecies));
%     stem3(temp1,temp2,temp3,'filled');
    
%     % plot momentum if required
%     if(strcmp(type,'rv'))
% 
%         %plot data for v
%         cutoff=optsPlot.vCutoff *max(max(nRS));
%         mask=(nRS<cutoff);
%         meanVS1=meanVS(:,:,1);
%         meanVS2=meanVS(:,:,2);
% 
%         meanVS1(mask)=NaN;
%         meanVS2(mask)=NaN;
% 
%         hp1(iSpecies)=surf(hPa1(iSpecies),xRS(:,:,1),xRS(:,:,2),meanVS1);
%         hp2(iSpecies)=surf(hPa2(iSpecies),xRS(:,:,1),xRS(:,:,2),meanVS2);
%             
%         hold(hPa,'on')
% 
%     end

    % set output handles of surfaces
    handles.hr=hr;
    
%     handles.hp1=hp1;
%     handles.hp2=hp2;

    % set output handles of axes if not given
    if(isempty(plotHandles))
        handles.hRa=hRa;
        handles.hPa=hPa;
    end
        
end
