function convStruct = FexMatrices_Meanfield_Spherical(optsPhys,IDC)

    if(isfield(optsPhys,'nParticlesS'))
        nSpecies=length(optsPhys.nParticlesS);
    else
        nSpecies=optsPhys.nSpecies;
    end

    params = optsPhys.V2;
    
    f      = str2func(params.V2DV2);
    
    params  = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams = length(paramNames);

    
    if(isfield(optsPhys,'FexNum'))
        convTemp = IDC.ComputeConvolutionMatrix(@PhiUniversal,optsPhys.FexNum);  
    else
        convTemp = IDC.ComputeConvolutionMatrix(@PhiUniversal);  
    end

    for iS = 1:nSpecies
        for jS = 1:nSpecies
            convStruct(iS,jS).Conv = convTemp(:,:,iS,jS);  %#ok
        end
    end
       
%--------------------------------------------------------------------------
	function z = PhiUniversal(x,y) 
        % note that x should be a scalar and y a vector
        z = zeros(length(y),nSpecies,nSpecies);

        for iSpecies = 1:nSpecies
            for jSpecies = 1:nSpecies
                optsPhysIJ = getIJParams(iSpecies,jSpecies);

                z(:,iSpecies,jSpecies) = f(x,y,optsPhysIJ);                
            end
        end                 
    end


    function optsPhysIJ = getIJParams(iS,jS)

        optsPhysIJ = optsPhys;
        
        for iParam=1:nParams
            
            paramValues = params.(paramNames{iParam});
            paramIJ = paramValues(iS,jS);
            optsPhysIJ.(paramNames{iParam}) = paramIJ;

        end           
    end

end
