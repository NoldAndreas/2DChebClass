function radialDistribution1D(X)
    
    Rij = getRij(X,X,1);
    
    Rij = triu(Rij);
    
    Rij = Rij(:);
    
    Rij = Rij(Rij>0);
    
    histogram(Rij)

end