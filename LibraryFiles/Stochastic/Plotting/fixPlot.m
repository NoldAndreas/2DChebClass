function fixPlot(h,xMin,xMax,yMin,yMax,xLab,yLab,time,legPos,legText)
% fixPlot(h,xMin,xMax,yMin,yMax,xLab,yLab,time,legPos,legText)
%   sets x and y axis limits and labels, adds a legend and prints the time
%   in the top right
%
% INPUTS:
%     h           -- axis handle
%     xMin        -- x axis minimum value
%     xMax        -- x axis maximum value
%     yMin        -- y axis minimum value
%     yMax        -- y axis maximum value
%     xLab        -- string for x axis label
%     yLab        -- string for y axis label
%     time        -- time to print
%     legPos      -- legend position, e.g. 'SouthEast'
%     legText     -- legend text in the form {'label1' 'label2'}

% set x and y limits of axes
xlim(h,[xMin, xMax]);
ylim(h,[yMin, yMax]);

% label axes
isLatex=~isempty(strfind(xLab,'$'));
if(isLatex)
    set(get(h,'XLabel'),'String',xLab,'Interpreter','LaTex','FontSize',20)
else
    set(get(h,'XLabel'),'String',xLab,'FontSize',20)
end

isLatex=~isempty(strfind(yLab,'$'));
if(isLatex)
    set(get(h,'YLabel'),'String',yLab,'Interpreter','LaTex','FontSize',20)
    set(get(h,'YLabel'),'Rotation',0.0)
else
    set(get(h,'YLabel'),'String',yLab,'FontSize',20)
end

% add a legend unless it's turned off

if(~strcmp(legPos,'off'))
    % determine if there is any latex in the string
    isLatex=~isempty(strfind(legText,'$'));
    % note that latex interpreter is slow and makes poor-looking text
    if(isLatex)
        legend(h,legText,'Location',legPos,'Interpreter','LaTex');
    else
        legend(h,legText,'Location',legPos);
    end
end

if(~isempty(time))
    % add the time to the plot
    axes(h);  %#ok % need this as you can't use gca with text
    textX=xMin+0.75*(xMax-xMin);
    textY=yMin+0.85*(yMax-yMin);
    text(textX,textY,['t=' num2str(time, '%6.4f')]);
end
    
% make sure that hold is off
hold(h,'off');
