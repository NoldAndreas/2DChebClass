function addEqLines(handlesRP,handlesM,equilibria,optsPlot)
%addEqLines(handlesRP,handlesM,equilibria)
%   add lines indicating equilibrium quantities to plots
%
% INPUTS: 
%  handlesRP  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  handlesM  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  equilbria -- structure of size [1 1] containing
%             REq             (vector of rho at positions xR)
%             vEq             (vector of v at positions xR)
%             xR              (vector positions, centres of histogram bins)
%             rMeanEq         (vector of mean positions, 
%                               same length as plotTimes)
%             pMeanEq         (vector of mean momenta, 
%                               same length as plotTimes)
%             plotTimes       (vector of plot times)

% if there are rho/v plots to do
if(~isempty(handlesRP))
    % get handles
    hRa=handlesRP.hRa;
    hPa=handlesRP.hPa;
    
    % get rho, v and histogram points
    REq=equilibria.REq;
    vEq=equilibria.vEq;
    xEq=equilibria.xEq;
    
    lineStyle=optsPlot.eqStyle{1};
    lineMarker=optsPlot.eqMarker{1};
    lineColour=optsPlot.eqColour{1};
    
    nSpecies=size(REq,2);
    for iSpecies=1:nSpecies
        hold(hRa,'on')
        
        if(optsPlot.plotDensity && strcmp(optsPlot.geom,'spherical'))
            nRdens=REq(:,iSpecies)./(4*pi*xEq(:,iSpecies).^2);
            hr=plot(hRa,xEq(:,iSpecies),nRdens);
        else
            
            hr=plot(hRa,xEq(:,iSpecies),REq(:,iSpecies));
        end
        set(hr,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
        
        hold(hPa,'on')
        hp=plot(hPa,xEq(:,iSpecies),vEq(:,iSpecies));
        set(hp,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
    end
end

% if there are mean r/p plots to do
if(~isempty(handlesM))
    % get handles
    hMRa=handlesM.hRa;
    hMPa=handlesM.hPa;

    % get means and times
    rMeanEq=equilibria.rMeanEq;
    vMeanEq=equilibria.vMeanEq;
    plotTimes=equilibria.plotTimes;

    handles.hRa=hMRa;
    handles.hPa=hMPa;
    plotLimit=length(plotTimes);
    type='rv';
    
    % eqStyle is of the form {{':','-'}} to be the same as the other styles
    optsPlot.lineStyle=optsPlot.eqStyle{1};
    optsPlot.lineMarker=optsPlot.eqMarker{1};
    optsPlot.lineColour=optsPlot.eqColour{1};
    
    hold(hMRa,'on');
    hold(hMPa,'on');
    plotRVmeans(rMeanEq,vMeanEq,plotTimes,plotLimit,optsPlot,handles,type);
    
    hold(hMRa,'off')
    hold(hMPa,'off')
end

