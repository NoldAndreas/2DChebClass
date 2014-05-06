function handles=plotRVmeans(meanr,meanv,plotTimes,plotLimit,optsPlot,handles,type)
%plotRVmeans(meanr,meanv,plotTimes,plotLimit,optsPlot,handles,type)
%   plots mean position and mean momentum at times in plotTimes up to 
%   plotLimit in axes given by plotHandles
%
% INPUTS: 
%  meanr        -- (nPlots,nSpecies) matrix of mean positions
%  meanv        -- (nPlots,nSpecies) matrix of mean velocities
%  plotTimes    -- (nPlots,1) vector of plotting times
%  plotLimit    -- positive integer <= nPlots to plot up to
%  optsPlot     -- a structure of size [1 1] containing:
%                        lineStyle  (linestyles for stochastic plots, 
%                                    should be  of the form ['b  ';'--b'])
%  plotHandles  -- a structure of size [1 1] containing:
%                        hRa        (handle for position axis)
%                        hPa        (handle for momentum axis)
%  type         -- either 'r' or 'rv' to determine if we do velocity plots
%
% OUTPUTS:
%  handles      -- a structure containing the handles of the plot lines,
%                  hr for the position and if type=rv also hp for velocity

% if no figure/axes then create them
if(nargin<6)
    figure;
    hRa=subplot(2,1,1);
    hPa=subplot(2,1,2);
else
    hRa=handles.hRa;
    hPa=handles.hPa;
end

% range over which to plot
range=1:plotLimit;

% plotting options
lineStyle=optsPlot.lineStyle;
lineMarker=optsPlot.lineMarker;
lineColour=optsPlot.lineColour;

nSpecies=size(meanr,2);

hr=zeros(nSpecies,1);

if(strcmp(type,'rv'))
    hp=hr;   %#ok
end

for iSpecies=1:nSpecies

    lineStyleS=lineStyle{iSpecies};
    lineMarkerS=lineMarker{iSpecies};
    lineColourS=lineColour{iSpecies};

    % plot position mean
    hr(iSpecies)=plot(hRa,plotTimes(range),meanr(range,iSpecies));
    set(hr(iSpecies),'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);

    handles.hr=hr;

    hold(hRa,'on');
    
    % plot momentum mean if we calculate it
    if(strcmp(type,'rv'))
        hp=plot(hPa,plotTimes(range),meanv(range,iSpecies));
        set(hp,'LineStyle',lineStyleS,'Color',lineColourS,'Marker',lineMarkerS);
        handles.hp=hp;
        hold(hPa,'on');
    
    end
    
end

%plot(hRa,plotTimes,nSpecies^(-1)*sum(meanr,2),':k');

hold(hRa,'off');
hold(hPa,'off');

end