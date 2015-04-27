function [eta,deta_drho] = ShearViscosity(rho,opts)
    O = ones(size(rho));    
    if(isfield(opts,'viscosity'))
        optsShearViscosity = opts.viscosity;
    else
        optsShearViscosity = opts;
    end
    
    if(isfield(optsShearViscosity,'etaC'))
        eta        = optsShearViscosity.etaC*O;
        deta_drho  = 0*O;
    elseif(isfield(optsShearViscosity,'etaL1'))
        eta        = optsShearViscosity.etaL1*rho;
        deta_drho  = optsShearViscosity.etaL1*O;
    elseif(isfield(optsShearViscosity,'etaLiq') && isfield(optsShearViscosity,'etaVap'))
        etaLiq     = optsShearViscosity.etaLiq;
        etaVap     = optsShearViscosity.etaVap;
        rhoLiq     = opts.rhoLiq_sat;
        rhoVap     = opts.rhoGas_sat;
        
        eta        = etaVap + (rho-rhoVap)/(rhoLiq-rhoVap)*(etaLiq-etaVap);
        deta_drho  = O/(rhoLiq-rhoVap)*(etaLiq-etaVap);        
        
    end
end