function fixPlot2Dsurf(h,opts)
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

%--------------------------------------------------------------------------
% Get options
%--------------------------------------------------------------------------

% axis limits
xMin=opts.xMin;
xMax=opts.xMax;
yMin=opts.yMin;
yMax=opts.yMax;
zMin=opts.zMin;
zMax=opts.zMax;

% axis labels
xLab=opts.xLab;
yLab=opts.yLab;
zLab=opts.zLab;

% time
time=opts.time;

%--------------------------------------------------------------------------
% Set x, y and z limits of axes
%--------------------------------------------------------------------------

xlim(h,[xMin, xMax]);
ylim(h,[yMin, yMax]);
zlim(h,[zMin, zMax]);

%--------------------------------------------------------------------------
% Set view point
%--------------------------------------------------------------------------

viewPoint=opts.viewPoint;
set(h,'View',viewPoint);

%--------------------------------------------------------------------------
% Set aspect ratio
%--------------------------------------------------------------------------

pbaspect(h,[ (xMax-xMin)   (yMax-yMin)   max( (xMax-xMin) , (yMax-yMin) ) ]);

%--------------------------------------------------------------------------
% Label axes
%--------------------------------------------------------------------------

xticks=get(h,'XTick');
yticks=get(h,'YTick');


nx=length(xticks);
ny=length(yticks);
mx=floor(nx/2)+1;
my=floor(ny/2)+1;

if(nx>2)
    xLabelPos=[xticks(mx)+(xticks(mx+1)-xticks(mx))/2, -0.1*(yMax-yMin)+yMin, zMin];
end

if(ny>2)
    sgny=sign(viewPoint(1));

    if (sgny>0)
        xVal=0.1*(xMax-xMin)+xMax;
    else
        xVal=-0.1*(xMax-xMin)+xMin;
    end

    yLabelPos=[xVal, yticks(my)+(yticks(my+1)-yticks(my))/2, zMin];
end

% FIX THIS FOR VIEW FROM OTHER SIDE?


isLatex=~isempty(strfind(xLab,'$'));
if(isLatex)
    set(get(h,'XLabel'),'String',xLab,'Interpreter','LaTex','FontSize',20)
else
    set(get(h,'XLabel'),'String',xLab,'FontSize',20)
end

if(nx>2)
    set(get(h,'XLabel'),'Position',xLabelPos);
end

isLatex=~isempty(strfind(yLab,'$'));
if(isLatex)
    set(get(h,'YLabel'),'String',yLab,'Interpreter','LaTex','FontSize',20)
else
    set(get(h,'YLabel'),'String',yLab,'FontSize',20)
end

if(ny>2)
    set(get(h,'YLabel'),'Position',yLabelPos);
end

isLatex=~isempty(strfind(zLab,'$'));
if(isLatex)
    set(get(h,'ZLabel'),'String',zLab,'Interpreter','LaTex','FontSize',20)
else
    set(get(h,'ZLabel'),'String',zLab,'FontSize',20)
end

%--------------------------------------------------------------------------
% Add the time to the plot
%--------------------------------------------------------------------------

if(~isempty(time))
    axes(h);  %#ok % need this as you can't use gca with text
    
    textX=xMax;
    textY=yMin+0.75*(yMax-yMin);
    textZ=zMin+0.9*(zMax-zMin);
    
    text(textX,textY,textZ,['t=' num2str(time, '%6.4f')]);
end

%--------------------------------------------------------------------------
% Clip surfaces
%--------------------------------------------------------------------------

% find each surface in the current axes
hSurf = findobj(gca,'Type','surf');
nSurf=length(hSurf);

for iSurf=1:nSurf
    
    % get x, y and z data for each surface
    xData=get(hSurf(iSurf),'xdata');
    yData=get(hSurf(iSurf),'ydata');
    zData=get(hSurf(iSurf),'zdata');
    
    % find all x and y points outside plotting area
    mask= xData<xMin | xData>xMax | yData<yMin | yData>yMax;       
    
    % remove all z data outside plotting area
    zData(mask)=NaN;
    set(hSurf(iSurf),'zdata',zData); 
end

% make sure that hold is off
hold(h,'off');
