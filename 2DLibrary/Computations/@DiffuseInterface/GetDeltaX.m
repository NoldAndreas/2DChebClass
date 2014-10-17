 function deltaX = GetDeltaX(this,phi,theta)
    y2Max      = this.optsNum.PhysArea.y2Max;
    pt.y2_kv   =  y2Max;
    fsolveOpts = optimset('Display','off');    

    phiBorder = phi(this.IC.Ind.top);
    xBorder   = this.IC.GetCartPts.y1_kv(this.IC.Ind.top);

    [~,j] = min(abs(phiBorder));
    x1_ig = xBorder(j);

    [y1Int,~,flag] = fsolve(@phiX1,x1_ig,fsolveOpts);
    if(flag < 1)
        cprintf('*r','GetDeltaX: Finding phi(y1,y_2max)=0: No solution found.\n');
        deltaX = NaN;
    else
        deltaX    = y1Int - y2Max/tan(theta);
        disp(['Delta x = ',num2str(deltaX)]);
    end  

    function z = phiX1(y1)
        pt.y1_kv = y1;
        IP       = this.IC.SubShapePtsCart(pt);
        z        = IP*phi;
    end   
end      