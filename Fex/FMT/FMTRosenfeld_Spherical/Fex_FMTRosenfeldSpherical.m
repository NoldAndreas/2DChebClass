function FMT = Fex_FMTRosenfeldSpherical(rho,IntMatrFex,kBT)
    % get number of species
    nSpecies=size(rho,2);

    n0=zeros(size(rho,1),1);
    n1=n0;
    n2=n0;
    n3=n0;
    n1r=n0;
    n2r=n0;
    
    % Di gives the integration matrix for n_i with r denoting a vector n
    % FMT for multiple species works by summing the n_i over all species
    % (see e.g. Roth review paper)
    for iSpecies=1:nSpecies    
        n0  = n0 + IntMatrFex(iSpecies).D0 * rho(:,iSpecies);
        n1  = n1 + IntMatrFex(iSpecies).D1 * rho(:,iSpecies);
        n2  = n2 + IntMatrFex(iSpecies).D2 * rho(:,iSpecies);
        n3  = n3 + IntMatrFex(iSpecies).D3 * rho(:,iSpecies);
        n1r = n1r + IntMatrFex(iSpecies).D1r * rho(:,iSpecies);     
        n2r = n2r + IntMatrFex(iSpecies).D2r * rho(:,iSpecies);
    end
    

    % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i \delta \rho
    % where the integration matrix for \delta n_i \delta \rho is given by
    % Di.  The non-vector ones are the same, but the vector ones need to be
    % slightly modified due to lack of symmetry under change of sign of the
    % argument

    FMT = zeros(size(rho));
    
    for iSpecies = 1:nSpecies
        h = - IntMatrFex(iSpecies).D0   * log(1-n3);
        h = h + IntMatrFex(iSpecies).D1 * (n2./(1-n3)) ;
        h = h + IntMatrFex(iSpecies).D2 * ( n1./(1-n3) + (n2.^2 - n2r.^2 )./(8*pi*(1-n3).^2));
        h = h + IntMatrFex(iSpecies).D3 * ( n0./(1-n3) + (n1.*n2 - n1r.*n2r)./((1-n3).^2)  + (n2.^3 -3* n2.*n2r.*n2r)./(12*pi*(1-n3).^3) );
        h = h + IntMatrFex(iSpecies).D1rF * (-n2r./(1-n3));  
        h = h + IntMatrFex(iSpecies).D2rF * (- n1r./(1-n3) - (n2.*n2r)./(4*pi*(1-n3).^2) );
        FMT(:,iSpecies)=h;
    end
    
    % rescale for temperature
    FMT = FMT*kBT;
    
end