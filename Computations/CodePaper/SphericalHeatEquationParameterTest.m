close all

Nvals = 10:10:100;
Lvals = 1:2:20;

maxRelErr = zeros(length(Nvals),length(Lvals));
t         = maxRelErr;
DErr      = maxRelErr;
D2Err     = maxRelErr;
DDErr     = maxRelErr;
IntErr    = maxRelErr;

for iN = 1:length(Nvals)
    for iL = 1:length(Lvals)
        [maxRelErr(iN,iL),t(iN,iL),DErr(iN,iL),D2Err(iN,iL),DDErr(iN,iL),IntErr(iN,iL)] ...
            = SphericalHeatEquation(Nvals(iN),Lvals(iL));
    end
end


figure

imagesc(log10(maxRelErr))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);

textStrings = num2str(t(:),'%0.2f');  
textStrings = strtrim(cellstr(textStrings));  
[x,y] = meshgrid(1:length(Nvals),1:length(Lvals)); 
hStrings = text(x(:),y(:),textStrings(:),... 
                'HorizontalAlignment','center');
textColors = repmat(t(:) > 0,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(hStrings,{'FontSize'},num2cell(20*ones(size(hStrings)),2));

% figure
% 
% imagesc(log10(DErr))
% 
% colorbar
% 
% set(gca,'YTick',1:length(Nvals));
% set(gca,'YTickLabel',Nvals);
% set(gca,'XTick',1:length(Lvals));
% set(gca,'XTickLabel',Lvals);
% 
% figure
% 
% imagesc(log10(D2Err))
% 
% colorbar
% 
% set(gca,'YTick',1:length(Nvals));
% set(gca,'YTickLabel',Nvals);
% set(gca,'XTick',1:length(Lvals));
% set(gca,'XTickLabel',Lvals);
% 
% figure
% 
% imagesc(log10(DDErr))
% 
% colorbar
% 
% set(gca,'YTick',1:length(Nvals));
% set(gca,'YTickLabel',Nvals);
% set(gca,'XTick',1:length(Lvals));
% set(gca,'XTickLabel',Lvals);
% 
% figure
% 
% imagesc(log10(IntErr))
% 
% colorbar
% 
% set(gca,'YTick',1:length(Nvals));
% set(gca,'YTickLabel',Nvals);
% set(gca,'XTick',1:length(Lvals));
% set(gca,'XTickLabel',Lvals);
% 
% figure
% 
% imagesc(t)
% 
% colorbar
% 
% set(gca,'YTick',1:length(Nvals));
% set(gca,'YTickLabel',Nvals);
% set(gca,'XTick',1:length(Lvals));
% set(gca,'XTickLabel',Lvals);

