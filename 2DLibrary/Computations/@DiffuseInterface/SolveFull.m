function [rho,uv] = SolveFull(this,ic)

    M       = this.IC.M;
    Ind     = this.IC.Ind;
    Diff    = this.IC.Diff;
    Int     = this.IC.Int;
    y2Max   = this.optsNum.PhysArea.y2Max;
    PtsCart = this.IC.GetCartPts();
    
	IntBorderTop    = this.IC.borderTop.IntSc;
    IntBorderBottom = this.IC.borderBottom.IntSc;
    
    %IntNormal_Path   = this.IC.borderTop.IntNormal_Path;
    %ptsBorderTop     = this.IC.borderTop.Pts;
    OR                 = ones(sum(Ind.right),1);        
	%IntPathUpLow       = this.IC.borderTop.IntNormal + this.IC.borderBottom.IntNormal; %?? => to check!!
        
        
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
    
    Atau    = Cak*(eta*Diff.LapVec + (zeta+eta/3)*Diff.gradDiv);
    
    theta      = pi/2;
    nParticles = 0;
    y10        = 0;
    
    uwall                                    = [UWall*ones(M,1) ; zeros(M,1);];
    uvBound                                  = zeros(2*M,1);
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));                   

	x = fsolve(@f,ic);
    
    D_B = x(1);
    x   = x(2:end);
    rho = x([T;F;F]);
    uv  = x([F;T;T]);
    
    function y = f(x)
        
        D_B   = x(1);
        x     = x(2:end);
        
        rho_s = x([T;F;F]);
        uv_s  = x([F;T;T]);
                
        mu       = DoublewellPotential(rho_s,Cn) - Cn*Diff.Lap*rho_s;        
        
        uvBound([Ind.top;Ind.top]) = ...
               GetSeppecherSolutionCart([PtsCart.y1_kv(Ind.top)- y10,PtsCart.y2_kv(Ind.top)],UWall,D_A,D_B,theta);       
        
        rho_f = rho_s + rho_m;
        rho_f2 = [rho_f;rho_f];
        
        y   = [rho_f.*(Diff.div*uv_s) + (Diff.Dy1*rho_s).*uv_s([T;F]) + (Diff.Dy2*rho_s).*uv_s([F;T]);... %Diff.div*(rho_f2.*uv_s);...
               -rho_f2.*(Diff.grad*mu) + Atau*uv_s];   
           
        %BC at wall (y_2 = 0)
        y([Ind.bottom;F;F]) = Ind.normalBottom*(Diff.grad*rho_s) - g;
        y([F;Ind.bottom;F]) = uv_s([Ind.bottom;F]) - UWall;
        y([F;F;Ind.bottom]) = uv_s([Ind.bottom;F]);
        
        %BC at top (y_2 = y_{2,max})       
        %topDirection           = [cos(theta)*ETop,sin(theta)*ETop];
        %y([Ind.top;F;F])       = topDirection*(Diff.grad*rho_s); %TODO!!!
        
        Cmu            = - Diff.div*(rho_f2.*(Diff.grad*mu))...
                         + Cak*(zeta + 4/3*eta)*(Diff.Lap*Diff.div); 

    A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
    b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
    
        
        y([F;Ind.top;Ind.top]) = uv_s([Ind.top;Ind.top]) - uvBound([Ind.top;Ind.top]);
        
        %BC at left (y1 = -infinity)
        y([Ind.left;F;F]) = mu(Ind.left);
        y([F;Ind.left;F]) = uv_s([Ind.left;F]) - UWall;
        y([F;F;Ind.left]) = uv_s([Ind.left;F]);
        
        %BC at right (y1 = infinity)                        
        m12     = Cak*eta*[Diff.Dy2 , Diff.Dy1]*uv_s - Cn*(Diff.Dy1*rho_s).*(Diff.Dy2*rho_s);                                
        y([Ind.right;F;F]) = (OR*(IntBorderTop-IntBorderBottom))*m12 + ...
                            - mu(Ind.right).*(rho_s(Ind.right)+rho_m)*y2Max...
                            + mu(Ind.left).*(rho_s(Ind.left)+rho_m)*y2Max;        
        
        y([F;Ind.right;F]) = uv_s([Ind.right;F]) - UWall;
        y([F;F;Ind.right]) = uv_s([Ind.right;F]);      
        
        %Mass
        y   = [Int*rho_s - nParticles;y];
        %massFlux   = ((-1+rho_m)-(1+rho_m))*(y2Max-0)*UWall;      %due to mapping to infinity
        %u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,D_A,dB,theta);                        
        %rhoBorder  = InterpOntoBorder*(rho+rho_m);
        %m          = IntNormal_Path*(u_flow.*repmat(rhoBorder,2,1)) + massFlux;
        
    end

end