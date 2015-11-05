function handles=plotRhoVdistDDFT(rho,v,Interp,optsPlot,plotHandles,type)
%plotRhoVdistDDFT(rhov,r,optsPlot,plotHandles,type)
%   plots position and velocity distributions in axes given by plotHandles
%
% INPUTS: 
%  rhov         -- (2*N,1) vector with columns [rho; v] specified  
%                   at grid points r
%  r            -- (N,1) vector of grid points
%  optsPlot     -- a structure of size [1 1] containing:
%                        lineStyle  (linestyles for stochastic plots, 
%                                    should be  of the form ['b  ';'--b'])
%                        plotDensity (true/false whether to plot the
%                                        density; false->distruibution)
%  plotHandles  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  type         -- either 'r' or 'rv' determining if we do momentum plots
%
% OUTPUTS:
%  handles      -- a structure containing the handles of the plot lines,
%                  hr for the position and if type=rv also hp for velocity


% if no figure/axes, create them
if(isempty(plotHandles))
    figure;
    hRa=subplot(2,2,1);
    hPa=subplot(2,2,2);
else
    hRa=plotHandles.hRa;
    hPa=plotHandles.hPa;
end

nSpecies=size(rho,2);

% plotting options
lineStyle=optsPlot.lineStyle;
lineMarker=optsPlot.lineMarker;
lineColour=optsPlot.lineColour;

hr=zeros(nSpecies,1);
hp=hr;

% interpolation matrix
InterPol = Interp.InterPol;

% plotting points in physical space
y = Interp.pts;

for iSpecies=1:nSpecies

    lineStyleS=lineStyle{iSpecies};
    lineMarkerS=lineMarker{iSpecies};
    lineColourS=lineColour{iSpecies};
    
    rhoS=rho(:,iSpecies);
        
    % do interpolation
    rhoS = InterPol*rhoS;
   
    % plot distribution rather than density
    if(strcmp(optsPlot.geom,'spherical') && ~optsPlot.plotDensity)
        rhoS=rhoS.*y.^2*4*pi;       
    end
    
    hr(iSpecies)=plot(hRa,y,rhoS); 
    set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    hold(hRa,'on')
    
    % plot momentum if required
    if(strcmp(type,'rv'))
        % extract v from rhov
        vS=v(:,iSpecies);
        if(optsPlot.plotCurrent)
            % current j=rho*v
            if(strcmp(optsPlot.geom,'spherical'))
                % rho from DDFT is actually the density not the distribution
                rhoS=rhoS.*y.^2*4*pi;
            end
            jS=rhoS.*vS;
            hp(iSpecies)=plot(hPa,y,jS);            
        else
            %plot data for v
    %         % cut off small values of rho
    %         maxRho=max(rho);
    %         cutoff=optsPlot.vCutoff;
    %         mask=(rho>cutoff*maxRho);
    %         % plot data for (cut off) v
    %         hp=plot(hPa,r(mask),v(mask));

            vS = InterPol*vS;
    
            hp(iSpecies)=plot(hPa,y,vS);
        end
        
        set(hp(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);
        
        hold(hPa,'on')
        
    end

end

% set output
handles.hr=hr;
handles.hp=hp;