function y1Cart = ComputeInterfaceContourY2(this,level,y2)

	if(nargin < 3)
        y2       = this.y2;
    end
    if((nargin < 2) || isempty(level))
        level = 0.5;
    end
    
    f            = this.GetRhoEq();
    drho         = (this.optsPhys.rhoLiq_sat - this.optsPhys.rhoGas_sat);        
    rhoV         = this.optsPhys.rhoGas_sat + level*drho;
    
    fsolveOpts   = optimset('Display','off');
    y1Cart       = zeros(size(y2));
    
    y1I = 10;
    for i = 1:length(y2)
        pt.y2_kv  = y2(i);        
        
        %Get initial guess
        % (a) get collocation point with y2-value closest to y2(i)
        CartPts   = this.IDC.GetCartPts();
        [~,mark1] = min(abs(CartPts.y2_kv-y2(i)));
        mark1     = (CartPts.y2_kv == CartPts.y2_kv(mark1));
        % (b) get from all collocation points (y1,y2) the one which value
        % is closest to rhoV
        [~,mark2] = min(abs(f(mark1)-rhoV));
        y1I = CartPts.y1_kv(mark1);
        y1I = y1I(mark2);
        
        y1Cart(i) = fsolve(@rhoX1,y1I,fsolveOpts);
                
        if(this.optsNum.PhysArea.alpha_deg == 90)
            y1I = y1Cart(i);
        elseif(this.optsNum.PhysArea.alpha_deg < 90)
            y1I       = max(y1Cart(i),4);        
        else
            y1I       = min(y1Cart(i),-4); 
        end
    end
    
    if(nargin < 3)
        this.IsolineInterfaceY2 = y1Cart;
    end
    
    function z = rhoX1(y1)
        pt.y1_kv = y1;
        IP = this.IDC.SubShapePtsCart(pt);
        z  = IP*f-rhoV;
    end    

end