function addOneLeg(optsPlot,hFig)
%addOneLeg(optsPlot,hFig)
% Adds a single legend at the top of a figure.
%
% INPUTS:
%  optsPlot  -- a structure containing
%         legText       (legend text for each line, of the form 
%                           {label1,label2,...,labeln})
%         lineStyle     (line styles for each line, same form)              
%         lineMarker    (line markers for each line, same form)              
%         lineColour    (line colours for each line, same form)
%         perRow        (integer (max) number of labels per row of legend)
%  hFig     --  figure handle in which to draw legend

if(isfield(optsPlot,'lineColourStoc'))
    lineStyleStoc=optsPlot.lineStyleStoc;
    lineMarkerStoc=optsPlot.lineMarkerStoc;
    lineColourStoc=optsPlot.lineColourStoc;
else
    lineStyleStoc=[];
    lineMarkerStoc=[];
    lineColourStoc=[];
end

if(isfield(optsPlot,'lineColourDDFT'))
    lineStyleDDFT=optsPlot.lineStyleDDFT;
    lineMarkerDDFT=optsPlot.lineMarkerDDFT;
    lineColourDDFT=optsPlot.lineColourDDFT;
else
    lineStyleDDFT=[];
    lineMarkerDDFT=[];
    lineColourDDFT=[];
end

% get legend text, styles and colours
legText=optsPlot.legTextR;
lineStyleTemp=cat(2,lineStyleStoc,lineStyleDDFT);
lineMarkerTemp=cat(2,lineMarkerStoc,lineMarkerDDFT);
lineColourTemp=cat(2,lineColourStoc,lineColourDDFT);

% calculate the number of rows in the legend
nLines=size(legText,2);
perRow=optsPlot.perRow;
nRows=ceil(nLines/perRow);

% flatten the cell arrays of styles, markers and colours
% which are stored as {iCalc}{iSpecies}
nCalc=numel(lineStyleTemp);
lineStyleFull=[];
lineMarkerFull=[];
lineColourFull=[];

for iCalc=1:nCalc
    lineStyleFull=[lineStyleFull,lineStyleTemp{iCalc}];    %#ok
    lineMarkerFull=[lineMarkerFull,lineMarkerTemp{iCalc}]; %#ok
    lineColourFull=[lineColourFull,lineColourTemp{iCalc}]; %#ok
end

% switch to the desired figure
figure(hFig);

for iRow =1:nRows
    % create an invisible set of axes for each row; we'll use these to plot
    % the data and create a row legend
    hLa=axes('Visible','off','Position',[-1 -1 0.1 0.1],'HitTest','off');

    % find limits of line info for this row
    lineMin=(iRow-1)*perRow+1;
    lineMax=min(iRow*perRow,nLines);

    % plot a point in the invisible axes of the correct style and colour
    % for each line
    for iLine=lineMin:lineMax
        htemp=plot(hLa,0,0);         % arbitrary line
        lineStyle=lineStyleFull{iLine};
        lineMarker=lineMarkerFull{iLine};
        lineColour=lineColourFull{iLine};
        set(htemp,'LineStyle',lineStyle,'Color',lineColour,'Marker',lineMarker);
        hold on    
    end
    
    % create the legend for this row
    hL=legend(legText(lineMin:lineMax));

    % determine the offset for each legend (should be about the height of
    % the text)
    offset=0.025;
    
    % position it in the top centre of the figure.  Position is of the form
    % [x y width height] on a scale of [0 1].  If the width/height is too
    % small it's automatically centred about x,y so we deliberately choose
    % it too small.
    set(hL,'Position',[0.5,1-offset*iRow,0.01,0.01]);
    
    % change to a horizontal legend and turn the box off (or we get an
    % odd-sized box for each row legend)
    set(hL,'Orientation','Horizontal');
    set(hL,'Box','Off');
end