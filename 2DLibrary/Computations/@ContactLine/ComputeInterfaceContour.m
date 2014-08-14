function ComputeInterfaceContour(this,level)

    f            = this.rho_eq;    
    drho         = (this.optsPhys.rhoLiq_sat - this.optsPhys.rhoGas_sat);
    
    if((nargin < 2) || isempty(level))
        level = 0.5;
    end
    rhoV         = this.optsPhys.rhoGas_sat + level*drho;
    
    fsolveOpts   = optimset('Display','off');
    y2Cart       = zeros(size(this.y1));
    
    y2I = 4;
    for i = 1:length(this.y1)
        pt.y1_kv  = this.y1(i);        
        y2Cart(i) = fsolve(@rhoX2,y2I,fsolveOpts);                    
        y2I       = max(y2Cart(i),4);        
    end
    
    this.IsolineInterface = y2Cart;
    
    function z = rhoX2(y2)
        pt.y2_kv = y2;
        IP = this.HS.SubShapePtsCart(pt);
        z  = IP*f-rhoV;
    end    

end