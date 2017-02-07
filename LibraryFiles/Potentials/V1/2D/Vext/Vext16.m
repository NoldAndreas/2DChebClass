function [V,VDiff] = Vext16(y1,y2)
    %intBound is a structur with:
    %   - y1_l
    %   - y1_u
    %   - y2_l
    %   - y2_u

    d  = 0.2;
    
    [y2I,dy2Idy1,dy2ddy1] = Interface(y1);

    z           = (y2-y2I)/d;
    V           = tanh(z);      
    
    dVdy1       = -(1-(tanh(z)).^2).*dy2Idy1/d;    
    dVdy2       =  (1-(tanh(z)).^2)/d;    
    dVddy1      = -(1-(tanh(z)).^2).*dy2ddy1/d + ...
                   (-2*tanh(z).*(1-tanh(z).^2)).*(dy2Idy1/d).^2;     
    dVddy2      = -2*tanh(z).*(1-tanh(z).^2)/d^2;
    dVdy1dy2    = 2*tanh(z).*(1-tanh(z).^2).*dy2Idy1/d^2;
    
    VDiff      = struct('dy1',dVdy1,'dy2',dVdy2,...
                        'ddy1',dVddy1,'ddy2',dVddy2,'dy1dy2',dVdy1dy2);

    function [y2I,dy2Idy1,dy2ddy1] = Interface(y1)
        L       = 2;
        A       = 0.2;
        
        y2I     = A*sin(2*pi/L*y1);
        dy2Idy1 = A*cos(2*pi/L*y1)*2*pi/L;
        dy2ddy1 = -A*sin(2*pi/L*y1)*(2*pi/L)^2;
    end

end