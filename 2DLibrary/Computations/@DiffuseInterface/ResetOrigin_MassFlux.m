function rho = ResetOrigin_MassFlux(this,theta,rho)

    InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    ptsBorderTop     = this.IC.borderTop.Pts;
    
    UWall            = this.optsPhys.UWall;
    rho_m            = this.optsPhys.rho_m;
    y2Max            = this.optsNum.PhysArea.y2Max;
        
    ptsCartShift0    = this.IC.GetCartPts();
    %Solve for parameter D_B to ensure Mass Balance
    theta  = fsolve(@GetMassInflux,0);                 

    ptsCartShift0.y2_kv = ptsCartShift0.y2_kv + deltaY1;
    Interp              = this.IC.SubShapePtsCart(ptsCartShift0);  
	rho                 = Interp*rho;
    
    function m = GetMassInflux(theta)
            %Density at infinity : -1,
            %Density at -infinity: 1. 
            massFlux   = ((-1+rho_m)-(1+rho_m))*(y2Max-0)*UWall;      %due to mapping to infinity

            %********************************
            %previous version:
            %u_flow     = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
            %m          = IC.borderTop.IntNormal*(u_flow.*(repmat(rho,2,1)+rho_m)) + massFlux;
            %********************************

            u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,0,0,theta);                        
            
            ptsCartShift        = ptsCartShift0;
            ptsCartShift.y1_kv  = ptsCartShift.y1_kv - cos(;
            Interp              = this.IC.SubShapePtsCart(ptsCartShift);  
            rho_h               = Interp*rho;
                        
            rhoBorder  = InterpOntoBorder*(rho_h+rho_m);
            m          = IntNormal_Path*(u_flow.*repmat(rhoBorder,2,1)) + massFlux;
    end

end