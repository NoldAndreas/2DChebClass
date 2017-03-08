function handles=plotRhoVdistStoc(rho,v,boxes,optsPlot,plotHandles,type)
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

% plotting options
lineStyle=optsPlot.lineStyle;
lineMarker=optsPlot.lineMarker;
lineColour=optsPlot.lineColour;

nSpecies = length(optsPlot.nParticlesS);

% handles of lines
hr=zeros(nSpecies,1);
hp=hr;

for iSpecies=1:nSpecies
    
    rhoS=rho(:,iSpecies);
    vS=v(:,iSpecies);
    boxesS=boxes(:,iSpecies);
    
    lineStyleS=lineStyle{iSpecies};
    lineMarkerS=lineMarker{iSpecies};
    lineColourS=lineColour{iSpecies};
    
    if(strcmp(optsPlot.geom,'spherical'))
        % plot data for rho
        if(optsPlot.plotDensity)
            %fprintf('density\n');
            
            % compute and plot density
            rhoDens=rhoS./(4*pi*boxesS.^2);
            hr(iSpecies)=plot(hRa,boxesS,rhoDens);
        else
            % if we want to plot the distribution
            %fprintf('distribution\n');
            hr(iSpecies)=plot(hRa,boxesS,rhoS);
        end

        set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    else
        % just plot distribution   
        hr(iSpecies)=plot(hRa,boxesS,rhoS); 
        
        set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    end
    
    % set output
    handles.hr=hr;

    hold(hRa,'on')
    
    % plot momentum if required
    if(strcmp(type,'rv'))
        if(optsPlot.plotCurrent)
            % current j=rho*v
            jS=rhoS.*vS;
            hp=plot(hPa,boxesS,jS);
        else
            %plot data for v
            cutoff=optsPlot.vCutoff;
            mask=(rhoS>cutoff*max(rhoS));
            %hp(iSpecies)=plot(hPa,xRS(mask),meanPS(mask));
            hp(iSpecies)=plot(hPa,boxesS(mask),vS(mask));
        end

        set(hp(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);
        % set output
        handles.hp=hp;
        
        hold(hPa,'on')
    end

end
