function [y,phi,upsilon,J] = Fex_FMTRoth(x,IntMatrFex,kBT,RS)
 
    upsilon = []; % still need to think about this

    N    = size(IntMatrFex(1).AD.m0,2);
    N_AD = size(IntMatrFex(1).AD.m0,1);

    % get number of species
    nSpecies=size(x,2);

    if(nargout >= 4)
        fullOutput = true;  %this case is only implemented for one species!!
    else
        fullOutput = false;
    end
    
    if(size(x,1) > N)
        fullInput = true; % given rho, n1, n2, n1_v_1, n1_v_2
    else
        fullInput = false;
    end
    
    rho     = x(1:N,:);

    n0=zeros(size(IntMatrFex(1).AD.m0,1),1);
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
        R    = RS(iSpecies);
        
        if(fullInput)
            m0Add      = x((1:N_AD)+N,iSpecies);
            n2Add      = x((1:N_AD)+N+1*N_AD,iSpecies);
            m1_1Add    = x((1:N_AD)+N+2*N_AD,iSpecies);
            m1_2Add    = x((1:N_AD)+N+3*N_AD,iSpecies);
            m2_11Add   = x((1:N_AD)+N+4*N_AD,iSpecies);
            m2_22Add   = x((1:N_AD)+N+5*N_AD,iSpecies);
            m2_12Add   = x((1:N_AD)+N+6*N_AD,iSpecies);
        else
            m0Add       = IntMatrFex(iSpecies).AD.m0 * rhoS;
            n2Add       = IntMatrFex(iSpecies).AD.n2 * rhoS;
            m1_1Add     = IntMatrFex(iSpecies).AD.m1_1 * rhoS;
            m1_2Add     = IntMatrFex(iSpecies).AD.m1_2 * rhoS;
            m2_11Add    = IntMatrFex(iSpecies).AD.m2_11 * rhoS;
            m2_22Add    = IntMatrFex(iSpecies).AD.m2_22 * rhoS;
            m2_12Add    = IntMatrFex(iSpecies).AD.m2_12 * rhoS;
        end
        
        n0   = n0   +  m0Add/(2*pi*R) ;
        n2   = n2   +  n2Add ;
        
        m0   = m0   +  m0Add ;
        
        m1_1 = m1_1 + m1_1Add ;
        m1_2 = m1_2 + m1_2Add ;
         
        m2_11 = m2_11 + m2_11Add ;
        m2_22 = m2_22 + m2_22Add ;
        m2_12 = m2_12 + m2_12Add ;
        m2_21 = m2_21 + m2_12Add ;
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

    FMT = zeros(N,1);
    
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
        h =     IntMatrFex(iSpecies).AAD.m0 * ( phi.n0/(2*pi*R) + phi.m0 ) ;
        h = h + IntMatrFex(iSpecies).AAD.n2 * phi.n2 ;
        %h = h + IntMatrFex(iSpecies).AAD.m0 * phi.m0 ;
        
        h = h - IntMatrFex(iSpecies).AAD.m1_1 * phi.m1_1 ;  
        h = h - IntMatrFex(iSpecies).AAD.m1_2 * phi.m1_2 ;  
        
        h = h + IntMatrFex(iSpecies).AAD.m2_11 * phi.m2_11 ;  
        h = h + IntMatrFex(iSpecies).AAD.m2_22 * phi.m2_22 ;  
        h = h + IntMatrFex(iSpecies).AAD.m2_12 * ( phi.m2_12 + phi.m2_21 );  
        %h = h + IntMatrFex(iSpecies).AAD.m2_21 * phi.m2_21 ;  
        
        FMT(:,iSpecies)=h;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
    if(size(x,1) > N)
        % only implemented for one species!!
        y   = [FMT;...
               m0 - IntMatrFex(1).AD.m0 * rho;...
               n2 - IntMatrFex(1).AD.n2 * rho;...
               m1_1 - IntMatrFex(1).AD.m1_1 * rho;...
               m1_2 - IntMatrFex(1).AD.m1_2 * rho;...
               m2_11 - IntMatrFex(1).AD.m2_11 * rho;...
               m2_22 - IntMatrFex(1).AD.m2_22 * rho;...
               m2_12 - IntMatrFex(1).AD.m2_12 * rho;...
               ];
    else
        y = FMT;
    end
    
    if(~fullOutput)        
        J = [];
        return;
    end
    
    AAD_IM       = IntMatrFex(1).AAD;
    AD_IM        = IntMatrFex(1).AD;
    
    n2M1         = 1./(1-n2);
    n2M2         = 1./(1-n2).^2;
    n2M3         = 1./(1-n2).^3;
    
    J_FMT_Full_AADn0 = [zeros(N_AD),diag(n2M1),zeros(N_AD),zeros(N_AD),zeros(N_AD),zeros(N_AD),zeros(N_AD)];
    
    J_FMT_Full_AADm0 = [1/2/pi*C0 * diag(n2M1), 1/2/pi*C0 * diag(m0.*n2M2),zeros(N_AD),zeros(N_AD),zeros(N_AD),zeros(N_AD),zeros(N_AD)];
        
    J_FMT_Full_AADn2 = [ diag( n2M1/(2*pi*R) + C0*m0.*n2M2/(2*pi) ),...
                         diag( n0.*n2M2 + ( C0 * m0.^2 + C1 * m1_dot_m1 + C2 * m2_dot_m2 ).*n2M3/(2*pi) ),...
                         diag( C1 * m1_1.*n2M2/(2*pi) ),...
                         diag( C1 * m1_2.*n2M2/(2*pi) ),...
                         diag( C2 * m2_11.*n2M2/(2*pi) ),...
                         diag( C2 * m2_22.*n2M2/(2*pi) ),...
                         diag( C2 * m2_12.*n2M2/(2*pi) ) ];
    
    J_FMT_Full_AADm1_1 = [zeros(N_AD),...
                          diag( C1 * m1_1.*n2M2/(2*pi) ),...
                          diag( C1 * n2M1/(2*pi) ),...
                          zeros(N_AD), ...
                          zeros(N_AD), ...
                          zeros(N_AD), ...
                          zeros(N_AD)];
    
    J_FMT_Full_AADm1_2 = [zeros(N_AD),...
                          diag( C1 * m1_2.*n2M2/(2*pi) ),...
                          zeros(N_AD), ...
                          diag( C1 * n2M1/(2*pi) ),...
                          zeros(N_AD), ...
                          zeros(N_AD), ...
                          zeros(N_AD)];

    J_FMT_Full_AADm2_11 = [zeros(N_AD),...
                           diag( C2 * m2_11.*n2M2/(2*pi) ),...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           diag( C2 * n2M1/(2*pi) ),...
                           zeros(N_AD), ...
                           zeros(N_AD)];

    J_FMT_Full_AADm2_22 = [zeros(N_AD),...
                           diag( C2 * m2_22.*n2M2/(2*pi) ),...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           diag( C2 * n2M1/(2*pi) ),...
                           zeros(N_AD)];

    J_FMT_Full_AADm2_12 = [zeros(N_AD),...
                           diag( C2 * m2_12.*n2M2/(2*pi) ),...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           zeros(N_AD), ...
                           diag( C2 * n2M1/(2*pi))];

    J_FMT_Full_AADm2_21 = J_FMT_Full_AADm2_12;

    J_FMT_Full =   AAD_IM.m0 * (J_FMT_Full_AADn0/(2*pi*R)+...
                                    J_FMT_Full_AADm0)  ...
                 + AAD_IM.n2 * J_FMT_Full_AADn2  ...
                 - AAD_IM.m1_1 * J_FMT_Full_AADm1_1 ...
                 - AAD_IM.m1_2 * J_FMT_Full_AADm1_2 ...
                 + AAD_IM.m2_11 * J_FMT_Full_AADm2_11 ...
                 + AAD_IM.m2_22 * J_FMT_Full_AADm2_22 ...
                 + AAD_IM.m2_12 * ( J_FMT_Full_AADm2_12 + J_FMT_Full_AADm2_21 ); 
        
    J        = [zeros(N),J_FMT_Full*kBT;...
                -[AD_IM.m0;AD_IM.n2;AD_IM.m1_1;AD_IM.m1_2;AD_IM.m2_11;AD_IM.m2_22;AD_IM.m2_12], ...
                eye(7*N_AD)];
    
end