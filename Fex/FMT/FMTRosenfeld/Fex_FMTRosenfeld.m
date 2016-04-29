function [y,phi,upsilon,J] = Fex_FMTRosenfeld(x,IntMatrFex,kBT,RS)

    upsilon = []; % still need to think about this

    N    = size(IntMatrFex(1).AD.n1,2);
    N_AD = size(IntMatrFex(1).AD.n1,1);

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
    
    n1     = zeros(size(IntMatrFex(1).AD.n1,1),1);
    n2     = n1;
    n1_v_1 = n1;
    n1_v_2 = n1;
    n0     = n1; 
    
    % Di gives the integration matrix for n_i with r denoting a vector n
    % FMT for multiple species works by summing the n_i over all species
    % (see e.g. Roth review paper)
    for iSpecies=1:nSpecies    
        rhoS=rho(:,iSpecies);
        R    = RS(iSpecies);
        
        if(fullInput)
            n1Add      = x((1:N_AD)+N,iSpecies);
            n2Add      = x((1:N_AD)+N+N_AD,iSpecies);
            n1_v_1Add  = x((1:N_AD)+N+2*N_AD,iSpecies);
            n1_v_2Add  = x((1:N_AD)+N+3*N_AD,iSpecies);
        else
            n1Add       = IntMatrFex(iSpecies).AD.n1 * rhoS;
            n2Add       = IntMatrFex(iSpecies).AD.n2 * rhoS;
            n1_v_1Add   = IntMatrFex(iSpecies).AD.n1_v_1 * rhoS;
            n1_v_2Add   = IntMatrFex(iSpecies).AD.n1_v_2 * rhoS;    
        end
        
        n0   = n0   + n1Add/(2*pi*R);  % n0 = n1/(2*pi*R)
        n1   = n1   + n1Add;
        n2   = n2   + n2Add;
        
        n1_v_1 = n1_v_1Add;
        n1_v_2 = n1_v_2Add;         
    end
    
    % construct dot products
    n1_dot_n1 = n1_v_1.*n1_v_1 + n1_v_2.*n1_v_2;
    
    % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i / \delta \rho
    % where the integration matrix for \delta n_i / \delta \rho is given by
    % IntMatrFex.n_i.  The non-vector ones are the same, but the vector
    % ones need to have the opposite sign due to sign change replacing r-r'
    % with r'-r.

    FMT = zeros(N,1);
    
    phi.n0     = -log(1-n2);
    phi.n1     = 1/2/pi * n1./(1-n2);
    phi.n2     = n0./(1-n2) + 1/4/pi * (n1.^2 - n1_dot_n1 )./((1-n2).^2);
    phi.n1_v_1 = -1/2/pi * n1_v_1./(1-n2);
    phi.n1_v_2 = -1/2/pi * n1_v_2./(1-n2);

    
    for iSpecies = 1:nSpecies
        h =     IntMatrFex(iSpecies).AAD.n1 * (phi.n0/(2*pi*R) + phi.n1) ;
        %h = h + IntMatrFex(iSpecies).AAD.n1 * phi.n1 ;
        h = h + IntMatrFex(iSpecies).AAD.n2 * phi.n2 ;
        h = h - IntMatrFex(iSpecies).AAD.n1_v_1 * phi.n1_v_1 ;  
        h = h - IntMatrFex(iSpecies).AAD.n1_v_2 * phi.n1_v_2 ;  
        FMT(:,iSpecies)=h;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
    if(size(x,1) > N)
        % only implemented for one species!!
        y   = [FMT;...
               n1-IntMatrFex(1).AD.n1 * rho;...
               n2-IntMatrFex(1).AD.n2 * rho;...
               n1_v_1-IntMatrFex(1).AD.n1_v_1 * rho;...
               n1_v_2-IntMatrFex(1).AD.n1_v_2 * rho];
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

    
    J_FMT_Full_AADn0 = [zeros(N_AD),diag(n2M1),zeros(N_AD),zeros(N_AD)];
    
    
    J_FMT_Full_AADn1 = [1/2/pi * diag(n2M1), 1/2/pi * diag(n1.*n2M2),zeros(N_AD),zeros(N_AD)];    
    
    % n0 = n1/(2*pi*R)
    J_FMT_Full_AADn2 = [ diag( n2M1/(2*pi*R) + n1.*n2M2/(2*pi) ),...
                         diag( n0.*n2M2 + (n1.^2 - n1_dot_n1).*n2M3/(2*pi) ),...
                         diag( -n1_v_1.*n2M2/(2*pi) ),...
                         diag( -n1_v_2.*n2M2/(2*pi) ) ];
                     
    J_FMT_Full_AADn1_v_1 = [zeros(N_AD),...
                            diag( -n1_v_1.*n2M2/(2*pi) ),...
                            diag( -n2M1/(2*pi) ),...
                            zeros(N_AD)];
    
    J_FMT_Full_AADn1_v_2 = [zeros(N_AD),...
                            diag( -n1_v_2.*n2M2/(2*pi) ),...
                            zeros(N_AD), ...
                            diag( -n2M1/(2*pi) )
                            ];
    J_FMT_Full =   AAD_IM.n1 * (J_FMT_Full_AADn0/(2*pi*R)+...
                                    J_FMT_Full_AADn1)  ...
                 + AAD_IM.n2 * J_FMT_Full_AADn2  ...
                 - AAD_IM.n1_v_1 * J_FMT_Full_AADn1_v_1 ...
                 - AAD_IM.n1_v_2 * J_FMT_Full_AADn1_v_2;        
    
    J        = [zeros(N),J_FMT_Full*kBT;...
                -[AD_IM.n1;AD_IM.n2;AD_IM.n1_v_1;AD_IM.n1_v_2],eye(4*N_AD)];

    
    
end