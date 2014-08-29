 function deltaX = GetDeltaX(this,rho,theta)
    y2Max      = this.optsNum.PhysArea.y2Max;
    pt.y2_kv   =  y2Max;
    fsolveOpts = optimset('Display','off');    

    rhoBorder = rho(this.IC.Ind.top);
    xBorder   = this.IC.GetCartPts.y1_kv(this.IC.Ind.top);

    [~,j] = min(abs(rhoBorder));
    x1_ig = xBorder(j);

    [y1Int,~,flag] = fsolve(@rhoX1,x1_ig,fsolveOpts);
    if(flag < 1)
        cprintf('*r','GetDeltaX: Finding rho(y1,y_2max)=0: No solution found.\n');
        deltaX = NaN;
    else
        deltaX    = y1Int - y2Max*cos(theta);
        disp(['Delta x = ',num2str(deltaX)]);
    end  

    function z = rhoX1(y1)
        pt.y1_kv = y1;
        IP       = this.IC.SubShapePtsCart(pt);
        z        = IP*rho;
    end   
end      