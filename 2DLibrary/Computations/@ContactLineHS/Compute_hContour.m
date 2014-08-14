function Compute_hContour(this,level)

    rho_eq       = this.GetRhoEq(); 
    N            = this.y1_SpectralLine.N;
    drho         = (this.optsPhys.rhoLiq_sat - this.optsPhys.rhoGas_sat);
    
    if((nargin < 2) || isempty(level))
        level = 0.5;
    end
    rhoV         = this.optsPhys.rhoGas_sat + level*drho;
    
    fsolveOpts   = optimset('Display','off');
    y2Cart       = zeros(N,1);
    
    y2I = 10;
    for i = 1:N
        pt.y1_kv  = this.y1_SpectralLine.Pts.y(i);        
        y2Cart(i) = fsolve(@rhoX2,y2I,fsolveOpts);                    
        %y2I       = max(y2Cart(i),4);        
        y2I       = y2Cart(i);
    end
    
    this.hContour = y2Cart;
    
    function z = rhoX2(y2)
        pt.y2_kv = y2;
        IP = this.IDC.SubShapePtsCart(pt);
        z  = IP*rho_eq-rhoV;
    end    

end