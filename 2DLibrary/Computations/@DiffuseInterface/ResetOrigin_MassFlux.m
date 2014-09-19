function rho = ResetOrigin_MassFlux(this,theta,phi)

    InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    ptsBorderTop     = this.IC.borderTop.Pts;
    
    UWall            = this.optsPhys.UWall;
    phi_m            = this.optsPhys.phi_m;
    y2Max            = this.optsNum.PhysArea.y2Max;
        
    ptsCartShift0    = this.IC.GetCartPts();
    %Solve for parameter D_B to ensure Mass Balance
    theta  = fsolve(@GetMassInflux,0);                 

    ptsCartShift0.y2_kv = ptsCartShift0.y2_kv + deltaY1;
    Interp              = this.IC.SubShapePtsCart(ptsCartShift0);  
	phi                 = Interp*phi;
    
    function m = GetMassInflux(theta)
            %Density at infinity : -1,
            %Density at -infinity: 1. 
            massFlux   = ((-1+phi_m)-(1+phi_m))*(y2Max-0)*UWall;      %due to mapping to infinity

            %********************************
            %previous version:
            %u_flow     = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
            %m          = IC.borderTop.IntNormal*(u_flow.*(repmat(phi,2,1)+phi_m)) + massFlux;
            %********************************

            u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,0,0,theta);                        
            
            ptsCartShift        = ptsCartShift0;
            ptsCartShift.y1_kv  = ptsCartShift.y1_kv - cos(;
            Interp              = this.IC.SubShapePtsCart(ptsCartShift);  
            phi_h               = Interp*phi;
                        
            phiBorder  = InterpOntoBorder*(phi_h+phi_m);
            m          = IntNormal_Path*(u_flow.*repmat(phiBorder,2,1)) + massFlux;
    end

end