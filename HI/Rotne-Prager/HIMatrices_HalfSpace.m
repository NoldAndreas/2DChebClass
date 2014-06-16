function HIStruct = HIMatrices_HalfSpace(opts,IDC)

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
    
    weights = {'HIweight_11','HIweight_12','HIweight_21','HIweight_22'};
    
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

            % 1) make annulus
            shape.N    = paramsIJ.N;
            shape.L    = paramsIJ.L;
            shape.RMin = paramsIJ.RMin;
             
            area = InfAnnulus(shape);
            
            paramsIJ.HIfn = f11;
            
            HITemp11 = IDC.ComputeConvolutionFiniteSupport2(area,weights,IDC.Pts,paramsIJ);
            
            paramsIJ.HIfn = f12;
            
            HITemp12 = IDC.ComputeConvolutionFiniteSupport2(area,weights,IDC.Pts,paramsIJ);
            
            % padded with weight of 1, so start indexing at 2
            HIInt11 = [HITemp11(:,:,2), HITemp11(:,:,3) ; ...
                       HITemp11(:,:,4), HITemp11(:,:,5) ];

            HIInt12 = [HITemp12(:,:,2), HITemp12(:,:,3) ; ...
                       HITemp12(:,:,4), HITemp12(:,:,5) ];                   
                   
            HIStruct(iS,jS).HIInt11 = HIInt11;
            HIStruct(iS,jS).HIInt12 = HIInt12;
            
            HIStruct(jS,iS).HIInt11 = HIInt11;
            HIStruct(jS,iS).HIInt12 = HIInt12;
        end
    end

    %--------------------------------------------------------------------------

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