function PostProcess_FilmThickness(this)

%     if(nargin >1)
%         m       = 200;
%         y1      = y1Int(1) + (0:(m-1))'/(m-1)*(y1Int(2)-y1Int(1));
%         this.y1 = y1;
%         %Int_y1  = ([y1(2:end)-y1(1:end-1);0]/2 + [0;y1(2:end)-y1(1:end-1)])'/2;
%     else
        y1 = this.y1;
%    end

    rho      = this.rho_eq;    
    
	rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;     
    filmThickness  = zeros(size(y1));    
	y2Max          = this.optsNum.maxComp_y2 + 10;
    
    for iy1 = 1:length(y1)
        y1i = y1(iy1);
        filmThickness(iy1) = this.HS.doIntFLine([y1i y1i],[0.5 y2Max],rho-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
    end
    
    ComputeInterfaceContour(this,0.05);

    this.filmThickness  = filmThickness;
        
end