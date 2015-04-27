function [zeta,dzeta_drho] = BulkViscosity(rho,opts)
    O = ones(size(rho));
    
	if(isfield(opts,'viscosity'))
        optsBulkViscosity = opts.viscosity;
    else
        optsBulkViscosity = opts;
    end
    
    
    if(isfield(optsBulkViscosity,'zetaC'))
        zeta        = optsBulkViscosity.zetaC*O;
        dzeta_drho  = 0*O;
    end
end