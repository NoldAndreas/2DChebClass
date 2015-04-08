function [zeta,dzeta_drho] = BulkViscosity(rho,optsBulkViscosity)
    O = ones(size(rho));
    if(isfield(optsBulkViscosity,'zetaC'))
        zeta        = optsBulkViscosity.zetaC*O;
        dzeta_drho  = 0;
    end
end