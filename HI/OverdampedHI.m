function HI = OverdampedHI(rho,IntMatrHI,gradMu)
    % get number of species

    nSpecies=size(rho,2);
    N = size(rho,1);
    
    HI1=zeros(2*N,nSpecies);
    HI2=zeros(2*N,nSpecies);

    rhoRep = [rho;rho];
    
    rhoGradMu = rhoRep.*gradMu;
    
    for iSpecies=1:nSpecies
        for jSpecies=1:nSpecies
            W1 = IntMatrHI(iSpecies,jSpecies).W1;
            W2 = IntMatrHI(iSpecies,jSpecies).W2;
            
            HI1(:,iSpecies) = HI1(:,iSpecies) + W1*rhoRep(:,jSpecies);
            HI2(:,iSpecies) = HI2(:,iSpecies) + W2*rhoGradMu(:,jSpecies);
        end
    end    
    
    HI = HI1.*gradMu + HI2;
end
