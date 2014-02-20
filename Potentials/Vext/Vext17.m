function [V,VDiff] = Vext17(y1,y2)
        d  = 0.8;

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
            A       = 2;
            y2I0    = 4;
            
            y2I     = y2I0 + A*exp(-(y1/L).^2);
            dy2Idy1 = -2*y1/L^2*A.*exp(-(y1/L).^2);
            dy2ddy1 = (-2/L^2+(2*y1/L^2).^2)*A.*exp(-(y1/L).^2);
            
            dy2Idy1((y1 == inf) | (y1==-inf)) = 0;
            dy2ddy1((y1 == inf) | (y1==-inf)) = 0;
        end
end