function addEqLinesDDFT(handlesRP,handlesM,ddft,optsPlot,doI,doF)
%addEqLines(handlesRP,handlesM,equilibria)
%   add dotted lines indicating equilibrium quantities to plots
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
    
lineStyle=optsPlot.eqStyle{1};
lineMarker=optsPlot.eqMarker{1};
lineColour=optsPlot.eqColour{1};

ddft=ddft(1);

% get rho, r and interpolation matrices
rhoEqI=ddft.rhoI;
rhoEqF=ddft.rhoF;
y=ddft.shape.Pts.y;
y(abs(y)==inf)=0;
r=ddft.shape.Interp.pts;
InterPol=ddft.shape.Interp.InterPol;
Int=ddft.shape.Int;

plotTimes=optsPlot.plotTimes;

nSpecies=size(rhoEqI,2);

% if there are rho/v plots to do
if(~isempty(handlesRP))
    % get handles
    hRa=handlesRP.hRa;
    hPa=handlesRP.hPa;
            
    for iSpecies=1:nSpecies
    
        rhoEqIS=rhoEqI(:,iSpecies);
        rhoEqFS=rhoEqF(:,iSpecies);

        rhoEqIS=InterPol*rhoEqIS;
        rhoEqFS=InterPol*rhoEqFS;
             
        if(strcmp(optsPlot.geom,'spherical') && ~optsPlot.plotDensity)
            rhoEqIS=4*pi*rhoEqIS.*r.^2;
            rhoEqFS=4*pi*rhoEqFS.*r.^2;
        end
    
        vEqS=zeros(size(rhoEqIS));
    
        % plot rho
        hold(hRa,'on')
        if( doI )
            hI = plot(hRa,r,rhoEqIS);
            set(hI,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
        end
        if( doF )
            hF=plot(hRa,r,rhoEqFS);
            set(hF,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
        end
        
        % plot v
        hold(hPa,'on')
        hP=plot(hPa,r,vEqS);
        set(hP,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
    end
        
        
    hold(hRa,'off')
    hold(hPa,'off')
end

% if there are mean r/p plots to do
if(~isempty(handlesM))
    % get handles
    hMRa=handlesM.hRa;
    hMPa=handlesM.hPa;
          
    for iSpecies=1:nSpecies
  
        rhoEqIS=rhoEqI(:,iSpecies);
        rhoEqFS=rhoEqF(:,iSpecies);
        
        rEqIS = ( Int*(rhoEqIS.*y) ) ./ (Int*rhoEqIS);
        rEqFS = ( Int*(rhoEqFS.*y) ) ./ (Int*rhoEqFS);
        
        rEqIS=rEqIS*ones(size(plotTimes));
        rEqFS=rEqFS*ones(size(plotTimes));

        vEqS=zeros(size(rEqIS));

        % plot r means
        hold(hMRa,'on')
        if( doI )
            hI=plot(hMRa,plotTimes,rEqIS);
            set(hI,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
        end
        
        if( doF )
            hF=plot(hMRa,plotTimes,rEqFS);
            set(hF,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
        end
        
        % plot p means
        hold(hMPa,'on')
        hP=plot(hMPa,plotTimes,vEqS);
        set(hP,'LineStyle',lineStyle{iSpecies},'Color',lineColour{iSpecies},'Marker',lineMarker{iSpecies});
    end
    
    hold(hMPa,'off')
    hold(hMRa,'off')
    
end
