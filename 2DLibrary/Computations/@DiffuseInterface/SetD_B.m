function D_B = SetD_B(this,theta,rho,initialGuessDB)

    InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    ptsBorderTop     = this.IC.borderTop.Pts;
    D_A              = this.optsPhys.D_A;
    UWall            = this.optsPhys.UWall;
    rho_m            = this.optsPhys.rho_m;
    y2Max            = this.optsNum.PhysArea.y2Max;

    %Solve for parameter D_B to ensure Mass Balance
    D_B  = fsolve(@GetMassInflux,initialGuessDB);                 
    
    function m = GetMassInflux(dB)
            %Density at infinity : -1,
            %Density at -infinity: 1. 
            massFlux   = ((-1+rho_m)-(1+rho_m))*(y2Max-0)*UWall;      %due to mapping to infinity

            %********************************
            %previous version:
            %u_flow     = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
            %m          = IC.borderTop.IntNormal*(u_flow.*(repmat(rho,2,1)+rho_m)) + massFlux;
            %********************************

            u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,D_A,dB,theta);                        
            rhoBorder  = InterpOntoBorder*(rho+rho_m);

            m          = IntNormal_Path*(u_flow.*repmat(rhoBorder,2,1)) + massFlux;
    end

end