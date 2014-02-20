function y1Cart = ComputeInterfaceContourY2(this,level,y2)

	if(nargin < 3)
        y2           = this.y2;
    end
    if((nargin < 2) || isempty(level))
        level = 0.5;
    end
    
    f            = this.rho_eq;    
    drho         = (this.optsPhys.rhoLiq_sat - this.optsPhys.rhoGas_sat);
    
    
    rhoV         = this.optsPhys.rhoGas_sat + level*drho;
    
    fsolveOpts   = optimset('Display','off');
    y1Cart       = zeros(size(y2));
    
    y1I = 4;
    for i = 1:length(y2)
        pt.y2_kv  = y2(i);        
        y1Cart(i) = fsolve(@rhoX1,y1I,fsolveOpts);
        y1I       = max(y1Cart(i),4);        
    end
    
    if(nargin < 3)
        this.IsolineInterfaceY2 = y1Cart;
    end
    
    function z = rhoX1(y1)
        pt.y1_kv = y1;
        IP = this.HS.SubShapePtsCart(pt);
        z  = IP*f-rhoV;
    end    

end