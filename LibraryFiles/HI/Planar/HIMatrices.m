function HIStruct = HIMatrices(opts,IDC)

    optsPhys = opts.optsPhys;
    optsNum  = opts.optsNum;

    if(isfield(optsPhys,'nParticlesS'))
        nSpecies=length(optsPhys.nParticlesS);
    else
        nSpecies=optsPhys.nSpecies;
    end

    params = optsPhys.HI;
    optsNum = optsNum.HINum;
    
    if(isfield(optsNum,'HIPreprocess'))
        fPreprocess = str2func(optsNum.HIPreprocess);
        params = fPreprocess(params);
    end
    
    f11      = str2func(optsNum.HI11);
    f12      = str2func(optsNum.HI12);

    if(isfield(optsNum,'HIg'))
        g = str2func(optsNum.HIg);
    else
        g = str2func('g_id');
    end
    
    
    if(isfield(optsNum,'N'))
        params.N = optsNum.N;
    else
        params.N = IDC.N;
    end
    
    if(isfield(optsNum,'L'))
        params.L = optsNum.L;
    elseif(isfield(IDC,'L'))
        params.L = IDC.L;
    end
        
    paramNames = fieldnames(params);
    nParams = length(paramNames);
    
    HIStruct(nSpecies,nSpecies).HIInt11 = [];
    
    for iS = 1:nSpecies
        for jS = iS:nSpecies
            paramsIJ = getIJParams(iS,jS);
            HITemp11 = IDC.ComputeConvolutionMatrix(@F11,paramsIJ);  
            HITemp12 = IDC.ComputeConvolutionMatrix(@F12,paramsIJ);  
            HIStruct(iS,jS).HIInt11 = HITemp11;
            HIStruct(iS,jS).HIInt12 = HITemp12;
            
            HIStruct(jS,iS).HIInt11 = HITemp11;
            HIStruct(jS,iS).HIInt12 = HITemp12;
        end
    end

    %--------------------------------------------------------------------------
    function z = F11(r)         
        z = f11(r,paramsIJ).*g(r,paramsIJ);                    
    end

    function z = F12(r)         
        z = f12(r,paramsIJ).*g(r,paramsIJ);                    
    end

    function paramsIJ = getIJParams(iS,jS)

       paramsIJ = params;

        for iParam=1:nParams
            paramValues = params.(paramNames{iParam});
            
            if(size(paramValues,1)==nSpecies && size(paramValues,2)==nSpecies)
                paramIJ = paramValues(iS,jS,:);
                paramsIJ.(paramNames{iParam}) = paramIJ;
            end

        end           
    end

end