function Compute_DisjoiningPressure_IV(this)    
       
    y1        = this.y1_SpectralLine.Pts.y;
    epw       = this.optsPhys.V1.epsilon_w;
	rho_eq    = this.GetRhoEq();                
    p         = this.optsPhys.p;   
    kBT       = this.optsPhys.kBT;
    
    %Disjoining Pressure based on 2D density profile:    
    
    fB        = zeros(size(y1));            
    rho_ref   = kron(ones(this.IDC.Pts.N1,1),this.rho1D_wg);  
    
    
    Geometry = struct('R_in',0,'LR',2,...
                      'th1',0,'th2',pi,...
                      'N',[40,40],...
                      'Origin',[0,0.5]);    
        
    HS                = InfWedge(Geometry);
        
    %Phys_Area = struct('L1',2,'L2',2,'N',[51;50],'y2Min',0.5);
    PlotArea       = struct('y1Min',-5,'y1Max',5,'y2Min',0.5,'y2Max',5,...
                            'N2',80,'N1',80);

   % HS        = HalfSpace(Phys_Area);
    HS.InterpolationPlotCart(PlotArea,true);
    IntHS     = HS.ComputeIntegrationVector();
    
    HSPtsCart = HS.GetCartPts();
    val_BH2D  = epw*BarkerHenderson_2D(sqrt((HSPtsCart.y1_kv).^2 + (HSPtsCart.y2_kv+0.5).^2));    
    
    hw = waitbar(0,'Computing disjoining pressure IV');
    
    for iy1 = 1:length(y1)
        y1i = y1(iy1);
        
        pts       = HSPtsCart; 
        pts.y1_kv = pts.y1_kv + y1i;        
        IP        = this.IDC.SubShapePtsCart(pts);
        
        pts.y1_kv = y1i; 
        pts.y2_kv = 0.5;
        IP0       = this.IDC.SubShapePtsCart(pts);
        
%        fB(iy1) = IntHS*((IP*rho_eq).*val_BH2D) + kBT*(IP0*rho_eq) - p;
        fB(iy1) = IntHS*((IP*(rho_eq-rho_ref)).*val_BH2D) + kBT*(IP0*(rho_eq-rho_ref)) ;
        
        hw = waitbar(iy1/length(y1));
    end
    close(hw);
    
    this.disjoiningPressure_IV = fB;    
    
    SumRule_DisjoiningPressure(this,'IV');
end