function [eta,deta_drho] = ShearViscosity(rho,optsShearViscosity)
    O = ones(size(rho));    
    if(isfield(optsShearViscosity,'etaC'))
        eta        = optsShearViscosity.etaC*O;
        deta_drho  = 0;
    end
end