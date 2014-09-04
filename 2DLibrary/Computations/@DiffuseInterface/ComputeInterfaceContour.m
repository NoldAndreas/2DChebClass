 function y1Cart = ComputeInterfaceContour(this)
    
    y2           = this.IC.Pts.y2;
    rho          = this.rho; 

    fsolveOpts   = optimset('Display','off');
    y1Cart       = zeros(size(y2));

    y1I = 0;
    for i = 1:length(y2)
        pt.y2_kv  = y2(i);        
        y1Cart(i) = fsolve(@rhoX1,y1I,fsolveOpts);        
        y1I       = y1Cart(i);        
    end
    
	this.IsolineInterfaceY2 = y1Cart;    

    function z = rhoX1(y1)
        pt.y1_kv = y1;
        IP = this.IC.SubShapePtsCart(pt);
        z  = IP*rho;
    end   

end