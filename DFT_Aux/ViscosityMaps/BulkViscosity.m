function [zeta,dzeta_drho] = BulkViscosity(rho,opts)
    O = ones(size(rho));
    
	if(isfield(opts,'viscosity'))
        optsBulkViscosity = opts.viscosity;
    else
        optsBulkViscosity = opts;
    end
    
%     if(isfield(optsBulkViscosity,'zetaC'))
%         zeta        = optsBulkViscosity.zetaC*O;
%         dzeta_drho  = 0*O;
%     end   
        
    if(isfield(optsBulkViscosity,'zetaC'))
        zeta        = optsBulkViscosity.zetaC*O;
        dzeta_drho  = 0*O;
    elseif(isfield(optsBulkViscosity,'zetaL1'))
        zeta        = optsBulkViscosity.zetaL1*rho;
        dzeta_drho  = optsBulkViscosity.zetaL1*O;
    elseif(isfield(optsBulkViscosity,'zetaLiq') && isfield(optsBulkViscosity,'zetaVap'))
        zetaLiq     = optsBulkViscosity.zetaLiq;
        zetaVap     = optsBulkViscosity.zetaVap;
        rhoLiq     = opts.rhoLiq_sat;
        rhoVap     = opts.rhoGas_sat;
        
        zeta        = zetaVap + (rho-rhoVap)/(rhoLiq-rhoVap)*(zetaLiq-zetaVap);
        dzeta_drho  = O/(rhoLiq-rhoVap)*(zetaLiq-zetaVap);        
        
    end
    
end