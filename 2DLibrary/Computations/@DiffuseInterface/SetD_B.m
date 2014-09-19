function D_B = SetD_B(this,theta,phi,initialGuessDB)

    InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    ptsBorderTop     = this.IC.borderTop.Pts;
    D_A              = this.optsPhys.D_A;
    UWall            = this.optsPhys.UWall;
    phi_m            = this.optsPhys.phi_m;
    y2Max            = this.optsNum.PhysArea.y2Max;

    %Solve for parameter D_B to ensure Mass Balance
    D_B  = fsolve(@GetMassInflux,initialGuessDB);                 
    
    function m = GetMassInflux(dB)
            %Density at infinity : -1,
            %Density at -infinity: 1. 
            massFlux   = ((-1+phi_m)-(1+phi_m))*(y2Max-0)*UWall;      %due to mapping to infinity

            %********************************
            %previous version:
            %u_flow     = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
            %m          = IC.borderTop.IntNormal*(u_flow.*(repmat(phi,2,1)+phi_m)) + massFlux;
            %********************************

            u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,D_A,dB,theta);                        
            phiBorder  = InterpOntoBorder*(phi+phi_m);

            m          = IntNormal_Path*(u_flow.*repmat(phiBorder,2,1)) + massFlux;
    end

end