hf = figure('position',[0,0,1000,1000]);
ha = axes;

nDDFT = length(DDFTStruct);

for iDDFT = 1:nDDFT

    rhoI = DDFTStruct(iDDFT).dynamicsResult.rho_t(:,:,1);

    Interp = DDFTStruct(iDDFT).IDC.Interp.InterPol;

    rhoIInterp(:,:,iDDFT) = Interp*rhoI;
    
    N(iDDFT) = DDFTStruct(iDDFT).IDC.N1;
end

nSpecies = size(rhoIInterp,2);

colours = {'r','b','g'};
symbols = {'o','s','^'};


for iSpecies = 1:nSpecies
    L2norm = norm(rhoIInterp(:,iSpecies,end));
    for iDDFT = 2:nDDFT
        rhoIDiff = rhoIInterp(:,iSpecies,iDDFT) - rhoIInterp(:,iSpecies,iDDFT-1);
        L2Diff(iDDFT-1,iSpecies) = norm(rhoIDiff)/L2norm;
    end
                                     
    plot(N(2:end),L2Diff(:,iSpecies),['-',symbols{iSpecies},...
                                         colours{iSpecies}],...
                                         'MarkerSize',10,...
                                         'MarkerFaceColor',colours{iSpecies}); 
    hold on;
end

set(ha,'YScale','log');


%ylabel('$\frac{\| \rho_{N} - \rho_{N-2} \|_2}{\| \rho_{50} \|_2}$','interpreter','latex');
ylabel(['$\| \rho_{N} - \rho_{N-2} \|_2 \Big/ \| \rho_{' num2str(N(end)) '} \|_2$'],'interpreter','latex');
xlabel('$N$','interpreter','latex');

