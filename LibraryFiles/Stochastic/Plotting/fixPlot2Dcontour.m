function fixPlot2Dcontour(h,opts)

% axis limits
xMin=opts.xMin;
xMax=opts.xMax;
yMin=opts.yMin;
yMax=opts.yMax;

% axis labels
xLab=opts.xLab;
yLab=opts.yLab;

% time
time=opts.time;

% % legend
% legPos=opts.legPos;
% legText=opts.legText;

% view

% set x, y and z limits of axes
xlim(h,[xMin, xMax]);
ylim(h,[yMin, yMax]);

%--------------------------------------------------------------------------
% Set aspect ratio
%--------------------------------------------------------------------------

%pbaspect(h,[ (xMax-xMin)   (yMax-yMin)   1 ]);


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
    %set(get(h,'YLabel'),'Rotation',0.0)
else
    set(get(h,'YLabel'),'String',yLab,'FontSize',20)
end


% add a legend unless it's turned off

% if(~strcmp(legPos,'off'))
%     % determine if there is any latex in the string
%     isLatex=~isempty(strfind(legText,'$'));
%     % note that latex interpreter is slow and makes poor-looking text
%     if(isLatex)
%         legend(h,legText,'Location',legPos,'Interpreter','LaTex');
%     else
%         legend(h,legText,'Location',legPos);
%     end
% end

if(~isempty(time))
    % add the time to the plot
    axes(h);  %#ok % need this as you can't use gca with text
    
    textX=0.65*xMax;
    textY=yMin+0.95*(yMax-yMin);
    
    text(textX,textY,['t=' num2str(time, '%6.4f')]);
end
  
% make sure that hold is off
hold(h,'off');
