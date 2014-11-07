function Compute_DisjoiningPressure_II(this)    
       
    y1  = this.y1_SpectralLine.Pts.y;

    %Disjoining Pressure based on 2D density profile:    
    rho_eq    = this.GetRhoEq();    
    %rho_ref   = kron(ones(this.IDC.Pts.N1,1),this.rho1D_wl);
    rho_ref   = kron(ones(this.IDC.Pts.N1,1),this.rho1D_wg);       
    
    fB            = zeros(size(y1));           
    
    y2Max = this.optsNum.maxComp_y2 + 10;
    
    hw = waitbar(0,'Computing disjoining pressure II');
    for iy1 = 1:length(y1)
        y1i = y1(iy1);
        
        pts.y1_kv = y1i; 
        pts.y2_kv = 0.5;
        IP0       = this.IDC.SubShapePtsCart(pts);

        [I,w,weights,IP,pts] = this.IDC.doIntFLine([y1i y1i],[0.5 y2Max],[],'CHEB');
        [h1,dVadd_i]         = getVAdd(pts.y1_kv,pts.y2_kv,0,this.optsPhys.V1);
        
        fB(iy1)              = - weights*(dVadd_i.dy2.*(IP*(rho_eq-rho_ref))) + this.optsPhys.kBT*(IP0*(rho_eq - rho_ref));        
        
        hw = waitbar(iy1/length(y1));
    end
    close(hw);        
        
    this.disjoiningPressure_II = fB;    
    
    SumRule_DisjoiningPressure(this,'II');
            
end