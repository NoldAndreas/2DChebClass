function [FMT,phi,upsilon] = Fex_FMTRosenfeld_3DFluid(rho,IntMatrFex,kBT,RS)
    % get number of species
    nSpecies=size(rho,2);

    n0      = zeros(size(IntMatrFex(1).AD.n2,1),1);
    n1      = n0;
    n2      = n0;
    n3      = n0;
    n1_v_1  = n0;
    n1_v_2  = n0;
	n2_v_1  = n0;
    n2_v_2  = n0;
    
    % Di gives the integration matrix for n_i with r denoting a vector n
    % FMT for multiple species works by summing the n_i over all species
    % (see e.g. Roth review paper)
    
    for iSpecies=1:nSpecies    
        rhoS = rho(:,iSpecies);
        R    = RS(iSpecies);
        
        n2Add       = IntMatrFex(iSpecies).AD.n2 * rhoS;
        n2_v_1Add   = IntMatrFex(iSpecies).AD.n2_v_1 * rhoS;
        n2_v_2Add   = IntMatrFex(iSpecies).AD.n2_v_2 * rhoS;
        
        n0     = n0   +  n2Add/(4*pi*R^2);
        n1     = n1   +  n2Add/(4*pi*R);        
        n2     = n2   +  n2Add;
        n3     = n3   +  IntMatrFex(iSpecies).AD.n3 * rhoS ;
        
        n1_v_1 = n1_v_1 + n2_v_1Add/(4*pi*R);
        n1_v_2 = n1_v_2 + n2_v_2Add/(4*pi*R);
        
        n2_v_1 = n2_v_1 + n2_v_1Add;
        n2_v_2 = n2_v_2 + n2_v_2Add;
         
    end
    
    % construct dot products
    n2_dot_n2 = n2_v_1.*n2_v_1 + n2_v_2.*n2_v_2;
    n1_dot_n2 = n1_v_1.*n2_v_1 + n1_v_2.*n2_v_2;
    
    % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i / \delta \rho
    % where the integration matrix for \delta n_i / \delta \rho is given by
    % IntMatrFex.n_i.  The non-vector ones are the same, but the vector
    % ones need to have the opposite sign due to sign change replacing r-r'
    % with r'-r.

    FMT = zeros(size(rho));
    
    phi.n0     = -log(1-n3);
    phi.n1     = n2./(1-n3);
    phi.n2     = n1./(1-n3) + 1/(8*pi) * (n2.^2 - n2_dot_n2 )./((1-n3).^2);
    phi.n3     = n0./(1-n3) + (n1.*n2 - n1_dot_n2 )./((1-n3).^2) +...
                            1/(12*pi) * (n2.^3 - 3*n2.*n2_dot_n2 )./((1-n3).^3);
    phi.n1_v_1 = -n2_v_1./(1-n3);
    phi.n1_v_2 = -n2_v_2./(1-n3);
    
    phi.n2_v_1 = -n1_v_1./(1-n3) - 1/(4*pi)*n2.*n2_v_1./((1-n3).^2);
    phi.n2_v_2 = -n1_v_2./(1-n3) - 1/(4*pi)*n2.*n2_v_2./((1-n3).^2);            
    
    for iSpecies = 1:nSpecies
        h = IntMatrFex(iSpecies).AAD.n2 * ...
                    (phi.n0/(4*pi*R^2) + phi.n1/(4*pi*R) + phi.n2);         
        h = h + IntMatrFex(iSpecies).AAD.n3 * phi.n3 ;         
         
        h = h - IntMatrFex(iSpecies).AAD.n2_v_1 * (phi.n1_v_1/(4*pi*R) + phi.n2_v_1 );  
        h = h - IntMatrFex(iSpecies).AAD.n2_v_2 * (phi.n1_v_2/(4*pi*R) + phi.n2_v_2 );                
                         
        FMT(:,iSpecies) = h;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
    %Compute grand potential
    if(nargout == 3)
        %first, add scalar
        upsilon = - n0.*log(1-n3) + n1.*n2./(1-n3) + (n2.^3)./(24*pi*(1-n3).^2);
        upsilon = upsilon -  n1_dot_n2./(1-n3) - (n2.*n2_dot_n2)./(8*pi*(1-n3).^2);        
        upsilon = upsilon*kBT;
    end
    
end


%*************************************************************************
%Old working version for multiple species
%(no speed up)
%*************************************************************************

% function [FMT,phi,upsilon] = Fex_FMTRosenfeld_3DFluid(rho,IntMatrFex,kBT)
%     % get number of species
%     nSpecies=size(rho,2);
% 
%     n0=zeros(size(IntMatrFex(1).AD.n0,1),1);
%     n1=n0;
%     n2=n0;
%     n3=n0;
%     n1_v_1=n0;
%     n1_v_2=n0;
% 	n2_v_1=n0;
%     n2_v_2=n0;
%     
%     % Di gives the integration matrix for n_i with r denoting a vector n
%     % FMT for multiple species works by summing the n_i over all species
%     % (see e.g. Roth review paper)
%     
%     for iSpecies=1:nSpecies    
%         rhoS=rho(:,iSpecies);
%         
%         n0   = n0   +  IntMatrFex(iSpecies).AD.n2 * rhoS ;  
%         n1   = n1   +  IntMatrFex(iSpecies).AD.n1 * rhoS ;
%         n2   = n2   +  IntMatrFex(iSpecies).AD.n2 * rhoS ;
%         n3   = n3   +  IntMatrFex(iSpecies).AD.n3 * rhoS ;
%         
%         n1_v_1 = n1_v_1 + IntMatrFex(iSpecies).AD.n1_v_1 * rhoS ;
%         n1_v_2 = n1_v_2 + IntMatrFex(iSpecies).AD.n1_v_2 * rhoS ;
%         
%         n2_v_1 = n2_v_1 + IntMatrFex(iSpecies).AD.n2_v_1 * rhoS ;
%         n2_v_2 = n2_v_2 + IntMatrFex(iSpecies).AD.n2_v_2 * rhoS ;
%          
%     end
%     
%     % construct dot products
%     n2_dot_n2 = n2_v_1.*n2_v_1 + n2_v_2.*n2_v_2;
%     n1_dot_n2 = n1_v_1.*n2_v_1 + n1_v_2.*n2_v_2;
%     
%     % now calculate \sum_i \int ( \partial Phi \partial n_i \delta n_i / \delta \rho
%     % where the integration matrix for \delta n_i / \delta \rho is given by
%     % IntMatrFex.n_i.  The non-vector ones are the same, but the vector
%     % ones need to have the opposite sign due to sign change replacing r-r'
%     % with r'-r.
% 
%     FMT = zeros(size(rho));
%     
%     phi.n0     = -log(1-n3);
%     phi.n1     = n2./(1-n3);
%     phi.n2     = n1./(1-n3) + 1/(8*pi) * (n2.^2 - n2_dot_n2 )./((1-n3).^2);
%     phi.n3     = n0./(1-n3) + (n1.*n2 - n1_dot_n2 )./((1-n3).^2) +...
%                             1/(12*pi) * (n2.^3 - 3*n2.*n2_dot_n2 )./((1-n3).^3);
%     phi.n1_v_1 = -n2_v_1./(1-n3);
%     phi.n1_v_2 = -n2_v_2./(1-n3);
%     
%     phi.n2_v_1 = -n1_v_1./(1-n3) - 1/(4*pi)*n2.*n2_v_1./((1-n3).^2);
%     phi.n2_v_2 = -n1_v_2./(1-n3) - 1/(4*pi)*n2.*n2_v_2./((1-n3).^2);            
%     
%     for iSpecies = 1:nSpecies
%         h =     IntMatrFex(iSpecies).AAD.n0 * phi.n0 ;
%         h = h + IntMatrFex(iSpecies).AAD.n1 * phi.n1 ;
%         h = h + IntMatrFex(iSpecies).AAD.n2 * phi.n2 ;
%         h = h + IntMatrFex(iSpecies).AAD.n3 * phi.n3 ;
%         h = h - IntMatrFex(iSpecies).AAD.n1_v_1 * phi.n1_v_1 ;  
%         h = h - IntMatrFex(iSpecies).AAD.n1_v_2 * phi.n1_v_2 ;  
%         h = h - IntMatrFex(iSpecies).AAD.n2_v_1 * phi.n2_v_1 ;  
%         h = h - IntMatrFex(iSpecies).AAD.n2_v_2 * phi.n2_v_2 ;  
%         FMT(:,iSpecies)=h;
%     end
%     
%     % rescale for temperature
%     FMT = real(FMT)*kBT;
%     
%     %Compute grand potential
%     if(nargout == 3)
%         %first, add scalar
%         upsilon = - n0.*log(1-n3) + n1.*n2./(1-n3) + (n2.^3)./(24*pi*(1-n3).^2);
%         upsilon = upsilon -  n1_dot_n2./(1-n3) - (n2.*n2_dot_n2)./(8*pi*(1-n3).^2);        
%         upsilon = upsilon*kBT;
%     end
%     
% end