function Compute_DisjoiningPressure_IV(this)    
       
    y1      = this.y1_SpectralLine.Pts.y;
    PtsCart = this.IDC.GetCartPts;
    Int     = this.IDC.Int;
    epw     = this.optsPhys.V1.epsilon_w;

    %Disjoining Pressure based on 2D density profile:    
    rho_eq    = this.GetRhoEq();            
    
    fB            = zeros(size(y1));                   
    
    
    Phys_Area = struct('L1',3,'L2',2,'N',[20;20],...
                       'y2Min',0);
    HS = HalfSpace(Phys_Area);
    
    hw = waitbar(0,'Computing disjoining pressure IV');
    
    
    for iy1 = 1:length(y1)
        y1i = y1(iy1);
        
        BH2Dx = epw*BarkerHenderson_2D(sqrt((PtsCart.y1_kv - y1i).^2 + (PtsCart.y2_kv+0.5).^2)); 
        
        fB(iy1) = Int*(rho_eq.*BH2Dx); 
        
        hw = waitbar(iy1/length(y1));
    end
    close(hw);        
        
    this.disjoiningPressure_IV = fB;    
    
    %TODO: subtract reference pressure, and check integrability of rho_eq.*BH2Dx
    
   % SumRule_DisjoiningPressure_II(this);
end