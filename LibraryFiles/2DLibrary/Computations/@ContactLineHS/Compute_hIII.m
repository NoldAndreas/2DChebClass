function hIII = Compute_hIII(this,rho)

	y1             = this.y1_SpectralLine.Pts.y;    
    
    if(nargin < 2)
        rho            = this.GetRhoEq();
    end
    
    if(this.IDC.alpha > pi/2)
        rho_ref = rho(this.IDC.Ind.top&this.IDC.Ind.right);
    else
        rho_ref = rho(this.IDC.Ind.top&this.IDC.Ind.left);
    end
    
	rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;     
    
    hIII           = zeros(size(y1));    
	y2Max          = this.optsNum.maxComp_y2 + 10;
    
    hw = waitbar(0,'Computing h_{III}');
    for iy1 = 1:length(y1)
        y1i       = y1(iy1);
        hIII(iy1) = abs(this.IDC.doIntFLine([y1i y1i],[0.5 y2Max],rho-rho_ref,'CHEB')/(rhoLiq_sat-rhoGas_sat));
        hw = waitbar(iy1/length(y1));
    end
    close(hw);

    if(nargin < 2)
        this.hIII  = hIII;   
    end
end