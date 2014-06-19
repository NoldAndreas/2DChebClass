function [rho,uv] = SolveFull(this,nParticles)

    M       = this.IC.M;
    Ind     = this.IC.Ind;
    Diff    = this.IC.Diff;
    Int     = this.IC.Int;
    y2Max   = this.optsNum.PhysArea.y2Max;
    PtsCart = this.IC.GetCartPts();
    
	%InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    %IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    %ptsBorderTop     = this.IC.borderTop.Pts;
        
    D_A     = this.optsPhys.D_A;
    UWall   = this.optsPhys.UWall;
    Cn      = this.optsPhys.Cn;
    Cak     = this.optsPhys.Cak;
	eta     = this.optsPhys.eta;
    zeta    = this.optsPhys.zeta;
    rho_m   = this.optsPhys.rho_m;    	
    g       = this.optsPhys.g;    

    F       = false(M,1);
    T       = true(M,1);
    E       = eye(M);
	ETop    = E(Ind.top,:);            
    
    uwall                                    = [UWall*ones(M,1) ; zeros(M,1);];
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));                   

	x = fsolve(@f);
    
    function y = f(x)
        
        D_B   = x(1);
        x     = x(2:end);
        
        rho_s = x([T;F;F]);
        uv_s  = x([F;T;T]);
        
        mu    = DoublewellPotential(rho_s,Cn) - Cn*Diff.Lap*rho_s;
        theta = ...;         
        
        uvBound([Ind.top;Ind.top]) = ...
               GetSeppecherSolutionCart([PtsCart.y1_kv(Ind.top)- y10,PtsCart.y2_kv(Ind.top)],UWall,D_A,D_B,theta);       
        
        rho_f = rho_s + rho_m;
        
        y   = [rho_f*Diff.div*uv_s;...
               -rho_f*Diff.grad*mu/Cak + eta*Diff.LapVec*uv_s + (zeta+eta/3)*Diff.gradDiv*uv_s];                   
           
        %BC at wall (y_2 = 0)
        y([Ind.bottom;F;F]) = Ind.normalBottom*(Diff.grad*rho_s) - g;
        y([F;Ind.bottom;F]) = uv_s([F;Ind.bottom;F]) - UWall;
        y([F;F;Ind.bottom]) = uv_s([F;Ind.bottom;F]);
        
        %BC at top (y_2 = y_{2,max})
        
        topDirection           = [cos(theta)*ETop,sin(theta)*ETop];
        y([Ind.top;F;F])       = topDirection*(Diff.grad*rho_s); %Ind.normalTop*(Diff.grad*rho_s); %TODO!!!
        y([F;Ind.top;Ind.top]) = uv_s([F;Ind.top;Ind.top]) - uvBound([Ind.top;Ind.top]);
        
        %BC at left (y1 = -infinity)
        y([Ind.left;F;F]) = Diff.Dy2(Ind.left,:)*mu;
        y([F;Ind.left;F]) = uv_s([F;Ind.bottom;F]) - UWall;
        y([F;F;Ind.left]) = uv_s([F;Ind.bottom;F]);
        
        %BC at right (y1 = infinity)
        y([Ind.right;F;F]) = Diff.Dy2(Ind.right,:)*mu;
        y([F;Ind.right;F]) = uv_s([F;Ind.bottom;F]) - UWall;
        y([F;F;Ind.right]) = uv_s([F;Ind.bottom;F]);
        
        %Mass
        y   = [Int*rho_s - nParticles;y];
        %massFlux   = ((-1+rho_m)-(1+rho_m))*(y2Max-0)*UWall;      %due to mapping to infinity
        %u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,D_A,dB,theta);                        
        %rhoBorder  = InterpOntoBorder*(rho+rho_m);
        %m          = IntNormal_Path*(u_flow.*repmat(rhoBorder,2,1)) + massFlux;
        
    end

end