function FMT = Fex_FMTRosenfeld(rho,IntMatrFex,kBT,R)
    % get number of species
    nSpecies=size(rho,2);

    n0=zeros(size(IntMatrFex(1).AD.n0,1),1);
    n1=n0;
    n2=n0;
    n1_v_1=n0;
    n1_v_2=n0;
    
    % Di gives the integration matrix for n_i with r denoting a vector n
    % FMT for multiple species works by summing the n_i over all species
    % (see e.g. Roth review paper)
    for iSpecies=1:nSpecies    
        rhoS=rho(:,iSpecies);
        
        n0   = n0   +  IntMatrFex(iSpecies).AD.n0 * rhoS ;
        n1   = n1   +  IntMatrFex(iSpecies).AD.n1 * rhoS ;
        n2   = n2   +  IntMatrFex(iSpecies).AD.n2 * rhoS ;
        
        n1_v_1 = n1_v_1 + IntMatrFex(iSpecies).AD.n1_v_1 * rhoS ;
        n1_v_2 = n1_v_2 + IntMatrFex(iSpecies).AD.n1_v_2 * rhoS ;         
    end
    
    % construct dot products
    n1_dot_n1 = n1_v_1.*n1_v_1 + n1_v_2.*n1_v_2;
    
    % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i / \delta \rho
    % where the integration matrix for \delta n_i / \delta \rho is given by
    % IntMatrFex.n_i.  The non-vector ones are the same, but the vector
    % ones need to have the opposite sign due to sign change replacing r-r'
    % with r'-r.

    FMT = zeros(size(rho));
    
    phi.n0     = -log(1-n2);
    phi.n1     = 1/2/pi * n1./(1-n2);
    phi.n2     = n0./(1-n2) + 1/4/pi * (n1.^2 - n1_dot_n1 )./((1-n2).^2);
    phi.n1_v_1 = -1/2/pi * n1_v_1./(1-n2);
    phi.n1_v_2 = -1/2/pi * n1_v_2./(1-n2);

    
    for iSpecies = 1:nSpecies
        h =     IntMatrFex(iSpecies).AAD.n0 * phi.n0 ;
        h = h + IntMatrFex(iSpecies).AAD.n1 * phi.n1 ;
        h = h + IntMatrFex(iSpecies).AAD.n2 * phi.n2 ;
        h = h - IntMatrFex(iSpecies).AAD.n1_v_1 * phi.n1_v_1 ;  
        h = h - IntMatrFex(iSpecies).AAD.n1_v_2 * phi.n1_v_2 ;  
        FMT(:,iSpecies)=h;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
end