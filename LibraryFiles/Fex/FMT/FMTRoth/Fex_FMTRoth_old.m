function FMT = Fex_FMTRoth(rho,IntMatrFex,kBT,R)
    % get number of species
    nSpecies=size(rho,2);

    n0=zeros(size(IntMatrFex(1).AD.n0,1),1);
    n2=n0;
    m0=n0;
    
    m1_1=n0;
    m1_2=n0;
    
    m2_11=n0;
    m2_22=n0;
    m2_12=n0;
    m2_21=n0;
    
    % Di gives the integration matrix for n_i with r denoting a vector n
    % FMT for multiple species works by summing the n_i over all species
    % (see e.g. Roth review paper)
    for iSpecies=1:nSpecies    
        rhoS=rho(:,iSpecies);
        
        n0   = n0   +  IntMatrFex(iSpecies).AD.n0 * rhoS ;
        n2   = n2   +  IntMatrFex(iSpecies).AD.n2 * rhoS ;
        
        m0   = m0   +  IntMatrFex(iSpecies).AD.m0 * rhoS ;
        
        m1_1 = m1_1 + IntMatrFex(iSpecies).AD.m1_1 * rhoS ;
        m1_2 = m1_2 + IntMatrFex(iSpecies).AD.m1_2 * rhoS ;
         
        m2_11 = m2_11 + IntMatrFex(iSpecies).AD.m2_11 * rhoS ;
        m2_22 = m2_22 + IntMatrFex(iSpecies).AD.m2_22 * rhoS ;
        m2_12 = m2_12 + IntMatrFex(iSpecies).AD.m2_12 * rhoS ;
        m2_21 = m2_21 + IntMatrFex(iSpecies).AD.m2_21 * rhoS ;
    end
    
    % construct dot products
    m1_dot_m1 = m1_1.^2 + m1_2.^2; 
    m2_dot_m2 = m2_11.^2 + m2_22.^2 + m2_12.^2 + m2_21.^2;
    
    a = 11/4;
    %a = 1;  % to give Rosenfeld values
    
    C0 = (2+a)/3;
    C1 = (a-4)/3;
    C2 = (2-2*a)/3;
    
    % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i / \delta \rho
    % where the integration matrix for \delta n_i / \delta \rho is given by
    % IntMatrFex.n_i.  
    % The non-vector ones are the same, but the vector
    % ones need to have the opposite sign due to sign change replacing r-r'
    % with r'-r.
    % The matrix ones are also the same.

    FMT = zeros(size(rho));
    
    phi.n0     = -log(1-n2);
    phi.m0     = 1/2/pi * C0 * m0./(1-n2);
    phi.n2     = n0./(1-n2) + 1/4/pi * ( C0 * m0.^2 + C1 * m1_dot_m1 + C2 * m2_dot_m2 )./((1-n2).^2);
    
    phi.m1_1   = 1/2/pi * C1 * m1_1./(1-n2);
    phi.m1_2   = 1/2/pi * C1 * m1_2./(1-n2);
    
    phi.m2_11   = 1/2/pi * C2 * m2_11./(1-n2);
    phi.m2_22   = 1/2/pi * C2 * m2_22./(1-n2);
    phi.m2_12   = 1/2/pi * C2 * m2_12./(1-n2);
    phi.m2_21   = 1/2/pi * C2 * m2_21./(1-n2);
        
    for iSpecies = 1:nSpecies
        h =     IntMatrFex(iSpecies).AAD.n0 * phi.n0 ;
        h = h + IntMatrFex(iSpecies).AAD.n2 * phi.n2 ;
        h = h + IntMatrFex(iSpecies).AAD.m0 * phi.m0 ;
        
        h = h - IntMatrFex(iSpecies).AAD.m1_1 * phi.m1_1 ;  
        h = h - IntMatrFex(iSpecies).AAD.m1_2 * phi.m1_2 ;  
        
        h = h + IntMatrFex(iSpecies).AAD.m2_11 * phi.m2_11 ;  
        h = h + IntMatrFex(iSpecies).AAD.m2_22 * phi.m2_22 ;  
        h = h + IntMatrFex(iSpecies).AAD.m2_12 * phi.m2_12 ;  
        h = h + IntMatrFex(iSpecies).AAD.m2_21 * phi.m2_21 ;  
        
        FMT(:,iSpecies)=h;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
end