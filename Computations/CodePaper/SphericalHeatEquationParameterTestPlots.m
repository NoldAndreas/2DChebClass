close all
set(0,'defaultaxesfontsize',20);
set(0,'defaultlinelinewidth',2);

defaultPos = [0 0 1000 800];

saveDir = [pwd filesep 'Computations' filesep 'CodePaper' filesep 'Images' filesep];


hDErrDecay = figure('Position',defaultPos);
for iL=1:length(Lvals)
    set(gca,'LineStyleOrder','-|--|:');
    plot(Nvals,log10(DErr(:,iL)));
    hold all
end
legend(num2str(Lvals'))
xlabel('$N$','interpreter','latex');
ylabel('$\log ( \| Du - u''\| / \| u'' \| )$','interpreter','latex');

save2pdf([saveDir 'DErrDecay.pdf'],hDErrDecay);

hD2ErrDecay = figure('Position',defaultPos);
for iL=1:length(Lvals)
    set(gca,'LineStyleOrder','-|--|:');
    plot(Nvals,log10(D2Err(:,iL)));
    hold all
end
legend(num2str(Lvals'))
xlabel('$N$','interpreter','latex');
ylabel('$\log (\| D^2u - u''''\| / \| u'''' \| )$','interpreter','latex');
save2pdf([saveDir 'D2ErrDecay.pdf'],hD2ErrDecay);

hDDErrDecay = figure('Position',defaultPos);
for iL=1:length(Lvals)
    set(gca,'LineStyleOrder','-|--|:');
    plot(Nvals,log10(DDErr(:,iL)));
    hold all
end
legend(num2str(Lvals'))
xlabel('$N$','interpreter','latex');
ylabel('$\log (\| DDu - u''''\| / \| u'''' \|)$','interpreter','latex');
save2pdf([saveDir 'DDErrDecay.pdf'],hDDErrDecay);

hIntErrDecay = figure('Position',defaultPos);
for iL=1:length(Lvals)
    set(gca,'LineStyleOrder','-|--|:');
    plot(Nvals,log10(IntErr(:,iL) + 10^(-16)));
    hold all
end
legend(num2str(Lvals'))
xlabel('$N$','interpreter','latex');
ylabel('$\log (|Iu - \int u|/|\int u |)$','interpreter','latex');
save2pdf([saveDir 'IntErrDecay.pdf'],hIntErrDecay);



hmaxRelErr = figure('Position',defaultPos);

imagesc(log10(maxRelErr))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'maxRelErr.pdf'],hmaxRelErr);

% textStrings = num2str(t(:),'%0.2f');  
% textStrings = strtrim(cellstr(textStrings));  
% [x,y] = meshgrid(1:length(Lvals),1:length(Nvals)); 
% hStrings = text(x(:),y(:),textStrings(:),... 
%                 'HorizontalAlignment','center');
% textColors = repmat(ones(size(t(:))),1,3);  % white text
% set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
% 
% set(hStrings,{'FontSize'},num2cell(20*ones(size(hStrings)),2));



hDErr = figure('Position',defaultPos);

imagesc(log10(DErr))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'DErr.pdf'],hDErr);


hD2Err = figure('Position',defaultPos);

imagesc(log10(D2Err))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'D2Err.pdf'],hD2Err);



hDDErr=figure('Position',defaultPos);

imagesc(log10(DDErr))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'DDErr.pdf'],hDDErr);

hIntErr = figure('Position',defaultPos);

imagesc(log10(IntErr))

colorbar

set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'IntErr.pdf'],hIntErr);

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

