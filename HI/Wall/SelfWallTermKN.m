function [D,DD] = SelfWallTermKN(x,y,optsPhys)

    sigmaHS = optsPhys.sigmaHS;
    
    a = sigmaHS/2;
    
    yInv = y.^(-1);
    
    % no slip, keeping factor of D0 for main code
    % See Lauga & Squires, PoF, 17, 103102 (2005)
    
    Dx = (1 - 9/16*a*yInv + 1/8*(a*yInv).^3);
    Dy = (1 - 9/8*a*yInv + 1/2*(a*yInv).^3); % Note this becomes negative!!
    
    Dx_x = zeros(size(x));
    Dy_y = 9/8*a*(yInv).^2 - 3/2*a^3*yInv.^4;
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 

%     sigmaHS = optsPhys.sigmaHS;
%     
%     a = sigmaHS/2;
%     
%     yInv = y.^(-1);
%     
%     % From Kim and Netz, J Chem Phys 124 114709 (2006)
%     % Includes up to RP level which ensures positivity
%     
% %     Dx = (1 - 9/16*a*yInv + 1/8*(a*yInv).^3);
% %     Dy = (1 - 9/8*a*yInv + 1/2*(a*yInv).^3);
% %     
% %     Dx_x = zeros(size(x));
% %     Dy_y = 9/8*a*(yInv).^2 - 3/2*a^3*(yInv).^4;
% 
% 
%     Dx = (1-9/32*sigmaHS*yInv);
%     Dy = (1-9/16*sigmaHS*yInv); % Note this becomes negative!!
%     
%     Dx_x = zeros(size(x));
%     Dy_y = 9/16*sigmaHS*(yInv).^2;
% 
%     D  = [Dx;Dy];
%     DD = [Dx_x;Dy_y]; 

end