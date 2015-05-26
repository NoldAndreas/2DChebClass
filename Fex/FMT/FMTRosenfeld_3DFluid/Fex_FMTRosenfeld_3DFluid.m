function [y,phi,upsilon,J] = Fex_FMTRosenfeld_3DFluid(x,IntMatrFex,kBT,RS)%n_i,n2_i,n3_i,n2_v_1_i,n2_v_2_i

    N    = size(IntMatrFex(1).AD.n2,2); %size(n_i,1);
    N_AD = size(IntMatrFex(1).AD.n2,1); %size(n2_i,1);

    % get number of species
    nSpecies = size(x,2);    
    
    if(nargout >= 4)%size(x,1) > N)
        fullOutput = true;  %this case is only implemented for one species!!
    else
        fullOutput = false;
    end
    
    rho     = x(1:N,:);    
    %rhoS = n_i;
    
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
        
        if(fullOutput)
            n2Add      = x((1:N_AD)+N,iSpecies);
            n3Add      = x((1:N_AD)+N+N_AD,iSpecies);
            n2_v_1Add  = x((1:N_AD)+N+2*N_AD,iSpecies);
            n2_v_2Add  = x((1:N_AD)+N+3*N_AD,iSpecies);
        else
            n2Add       = IntMatrFex(iSpecies).AD.n2 * rhoS;
            n3Add       = IntMatrFex(iSpecies).AD.n3 * rhoS;
            
            if(isfield(IntMatrFex(iSpecies),'DiffAD'))
                n2_v_1Add   = -IntMatrFex(iSpecies).DiffAD.Dy1*n3Add;
                n2_v_2Add   = -IntMatrFex(iSpecies).DiffAD.Dy2*n3Add;
            else        
                n2_v_1Add   = IntMatrFex(iSpecies).AD.n2_v_1 * rhoS;
                n2_v_2Add   = IntMatrFex(iSpecies).AD.n2_v_2 * rhoS;
            end
                
        end                
        
        n0     = n0   +  n2Add/(4*pi*R^2);
        n1     = n1   +  n2Add/(4*pi*R);        
        n2     = n2   +  n2Add;
        n3     = n3   +  n3Add;
        
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

    %FMT = zeros(size(rho));
    FMT = zeros(N,1);
    
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
    
    %Compute grand potential
    if(nargout >= 3)
        %first, add scalar
        upsilon = - n0.*log(1-n3) + n1.*n2./(1-n3) + (n2.^3)./(24*pi*(1-n3).^2);
        upsilon = upsilon -  n1_dot_n2./(1-n3) - (n2.*n2_dot_n2)./(8*pi*(1-n3).^2);        
        upsilon = upsilon*kBT;
    end
    
    % rescale for temperature
    FMT = real(FMT)*kBT;
    
    if(size(x,1) > N)
        y   = [FMT;...
               n2-IntMatrFex(iSpecies).AD.n2 * rho;...
               n3-IntMatrFex(iSpecies).AD.n3 * rho;...
               n2_v_1-IntMatrFex(iSpecies).AD.n2_v_1 * rho;...
               n2_v_2-IntMatrFex(iSpecies).AD.n2_v_2 * rho];
    else
        y = FMT;
    end
    
    if(~fullOutput)        
        J = [];
        return;
    end
    
    %********************************
    % Compose full vector
    %********************************    
        
    %********************************
    % Compute Jacobian
    %********************************
    AAD_IM       = IntMatrFex(iSpecies).AAD;
    AD_IM        = IntMatrFex(iSpecies).AD;
    
    n3M1         = 1./(1-n3);
    n3M2         = 1./(1-n3).^2;
    n3M3         = 1./(1-n3).^3;
    n3M4         = 1./(1-n3).^4;
    
    J_FMT_Full_AADn0 = [zeros(N_AD),diag(n3M1),zeros(N_AD),zeros(N_AD)];
    J_FMT_Full_AADn1 = [diag(n3M1),diag(n2.*n3M2),zeros(N_AD),zeros(N_AD)];    
    J_FMT_Full_AADn2 = [diag(n3M1/(4*pi*R)+n2.*n3M2/(4*pi)),...
                        diag(n1.*n3M2 +(n2.^2 - n2_dot_n2).*n3M3/(4*pi)),...
                        diag(-n2_v_1.*n3M2/(4*pi)),...
                        diag(-n2_v_2.*n3M2/(4*pi))];
                    
    J_FMT_Full_AADn3 = [diag(n3M1/(4*pi*R^2) + n2.*n3M2/(4*pi*R) + n1.*n3M2 + (n2.^2 - n2_dot_n2).*n3M3/(4*pi)),...
                        diag(n0.*n3M2 + 2*(n1.*n2-n1_dot_n2).*n3M3+(n2.^3-3*n2.*n2_dot_n2).*n3M4/(4*pi)),...
                        diag( (-n2_v_1.*n3M2/(4*pi*R)-n1_v_1.*n3M2-n2.*n2_v_1.*n3M3/(2*pi))),...
                        diag( (-n2_v_2.*n3M2/(4*pi*R)-n1_v_2.*n3M2-n2.*n2_v_2.*n3M3/(2*pi)))];
                    
                    
    J_FMT_Full_AADn1_v_1 = [zeros(N_AD),...
                            diag(-n2_v_1.*n3M2),...
                            diag(-n3M1),...
                            zeros(N_AD)];
                        
    J_FMT_Full_AADn1_v_2 = [zeros(N_AD),...
                            diag(-n2_v_2.*n3M2),...                            
                            zeros(N_AD),...
                            diag(-n3M1)];                    
                    
    J_FMT_Full_AADn2_v_1 = [diag(-n2_v_1.*n3M2/(4*pi)),...
                            diag(-n1_v_1.*n3M2-n2.*n2_v_1.*n3M3/(2*pi)),...
                            diag(-n3M1/(4*pi*R)-n2.*n3M2/(4*pi)),...
                            zeros(N_AD)];
                        
    J_FMT_Full_AADn2_v_2 = [diag(-n2_v_2.*n3M2/(4*pi)),...
                            diag(-n1_v_2.*n3M2-n2.*n2_v_2.*n3M3/(2*pi)),...                            
                            zeros(N_AD),...
                            diag(-n3M1/(4*pi*R)-n2.*n3M2/(4*pi))];
                    
    J_FMT_Full =   AAD_IM.n2*(J_FMT_Full_AADn0/(4*pi*R^2)+...
                              J_FMT_Full_AADn1/(4*pi*R)  +...
                              J_FMT_Full_AADn2)  ...
                 + AAD_IM.n3*J_FMT_Full_AADn3  ...
                 - AAD_IM.n2_v_1*(J_FMT_Full_AADn1_v_1/(4*pi*R) + ...
                                  J_FMT_Full_AADn2_v_1) ...
                 - AAD_IM.n2_v_2*(J_FMT_Full_AADn1_v_2/(4*pi*R) + ...
                                  J_FMT_Full_AADn2_v_2);        
    
    J        = [zeros(N),J_FMT_Full;...
                -[AD_IM.n2;AD_IM.n3;AD_IM.n2_v_1;AD_IM.n2_v_2],eye(4*N_AD)];
            
    J = kBT*J;
        
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