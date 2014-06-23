function dFdRho=Fex_Percus(rho,Eta,kBT)

    nSpecies = size(rho,2);    
    rhoTotal = sum(rho,2);
    
    etaFrho = Eta.F * rhoTotal;    
    dFdRho = - log(1 - etaFrho) + Eta.B * ( rhoTotal .* (1 - etaFrho).^(-1) );
    dFdRho = kBT * dFdRho;
    
    dFdRho = repmat(dFdRho,1,nSpecies);

end