close all
set(0,'defaultaxesfontsize',20);
set(0,'defaultlinelinewidth',2);

defaultPos = [0 0 1000 800];

saveDir = [pwd filesep 'Computations' filesep 'AN_CodePaper' filesep 'Images' filesep];

Nvals = 10:10:200;
Lvals = 1:14;

IntErr = zeros(length(Nvals),length(Lvals));

for iN = 1:length(Nvals)
    iN
    for iL = 1:length(Lvals)
        IntErr(iN,iL) = IntegrationNL(Nvals(iN),Lvals(iL));
    end
end
        
hIntErr = figure('Position',defaultPos);
imagesc(log10(IntErr));
colormap(hIntErr,gray)
colorbar


set(gca,'YTick',1:length(Nvals));
set(gca,'YTickLabel',Nvals);
set(gca,'XTick',1:length(Lvals));
set(gca,'XTickLabel',Lvals);
xlabel('$L$','interpreter','latex');
ylabel('$N$','interpreter','latex');
save2pdf([saveDir 'IntErr.pdf'],hIntErr);

