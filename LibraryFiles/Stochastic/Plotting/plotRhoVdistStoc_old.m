function handles=plotRhoVdistStoc(x,p,optsPlot,plotHandles,type)
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

% if no figure/axes, create them
if(isempty(plotHandles))
    figure;
    hRa=subplot(2,1,1);
    hPa=subplot(2,1,2);
else
    hRa=plotHandles.hRa;
    hPa=plotHandles.hPa;
end

% determine the relevant parameters
geom=optsPlot.geom;
nBins=optsPlot.nBins;
dim=optsPlot.dim;

% get the R and P data
[R,P]=getRP(x,p,geom,dim);

nParticlesS=optsPlot.nParticlesS;
nSpecies=length(nParticlesS);

% % do histogramming of R and obtain v(r).  Note we can't just use hist as we want P=v
% % as a function of r, not as a histogram
[nR,meanP,xR]= RPhist(R,P,nBins,nParticlesS);

% convert to velocities
mS=optsPlot.mS;
mS=repmat(mS',size(meanP,1),1);
meanV=meanP./mS;

% plotting options
lineStyle=optsPlot.lineStyle;
lineMarker=optsPlot.lineMarker;
lineColour=optsPlot.lineColour;

% handles of lines
hr=zeros(nSpecies,1);
hp=hr;

for iSpecies=1:nSpecies
    
    nRS=nR(:,iSpecies);
    %meanPS=meanP(:,iSpecies);
    meanVS=meanV(:,iSpecies);
    xRS=xR(:,iSpecies);
    
    lineStyleS=lineStyle{iSpecies};
    lineMarkerS=lineMarker{iSpecies};
    lineColourS=lineColour{iSpecies};
    
    if(strcmp(optsPlot.geom,'spherical'))
        % plot data for rho
        if(optsPlot.plotDensity)
            %fprintf('density\n');
            
            % compute and plot density
            nRdens=nRS./(4*pi*xRS.^2);
            hr(iSpecies)=plot(hRa,xRS,nRdens);
        else
            % if we want to plot the distribution
            %fprintf('distribution\n');
            hr(iSpecies)=plot(hRa,xRS,nRS);
        end

        set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    else
        % just plot distribution   
        hr(iSpecies)=plot(hRa,xRS,nRS); 
        
%         hold(hRa,'on');
%         
%         D=1/2;
%         t0=-1/D;
%         
%         t=optsPlot.plotTime;
%         f=5*1/sqrt(pi*(t-t0)*D/0.1)*exp(-0.1*xRS.^2/D/(t-t0));
%         h=xRS(2)-xRS(1);
%         
%         sum(f)*h;
%         sum(nRS)*h;
%       
%         plot(hRa,xRS,f,'r'); 

        
        set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    end
    
    % set output
    handles.hr=hr;

    hold(hRa,'on')
    
    % plot momentum if required
    if(strcmp(type,'rv'))
        if(optsPlot.plotCurrent)
            % current j=rho*v
            jS=nRS.*meanVS;
            hp=plot(hPa,xRS,jS);
        else
            %plot data for v
            cutoff=optsPlot.vCutoff;
            mask=(nRS>cutoff*max(nRS));
            %hp(iSpecies)=plot(hPa,xRS(mask),meanPS(mask));
            hp(iSpecies)=plot(hPa,xRS(mask),meanVS(mask));
        end

        set(hp(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);
        % set output
        handles.hp=hp;
        
        hold(hPa,'on')
    end

end
