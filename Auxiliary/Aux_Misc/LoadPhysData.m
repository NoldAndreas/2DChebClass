    function [kBT,nParticles,HS_f] = LoadPhysData(optsPhys)        
        kBT        = optsPhys.kBT;   
        nParticles = optsPhys.nParticles;             
        if(isfield(optsPhys,'HSBulk'))
            HS_f       = str2func(optsPhys.HSBulk);
        end
    end
