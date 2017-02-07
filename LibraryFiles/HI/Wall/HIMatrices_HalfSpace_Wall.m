function HIStruct = HIMatrices_HalfSpace_Wall(opts,IDC)

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
    
    if(isfield(optsNum,'doConv'))
        doConv = optsNum.doConv;
    else
        doConv = false;
    end
    
    for iS = 1:nSpecies
        for jS = iS:nSpecies
            paramsIJ = getIJParams(iS,jS);

            % 1) make annulus
            shapeGeom.N     = paramsIJ.N;
            shapeGeom.L     = paramsIJ.L;
            shapeGeom.RMin  = paramsIJ.RMin;
            shapeGeom.shape = 'InfAnnulus';
            
            %HITemp11 = IDC.ComputeAreaIntegrationFiniteSupport(shapeGeom,@F11,[],doConv);
            %HITemp12 = IDC.ComputeAreaIntegrationFiniteSupport(shapeGeom,@F12,[],doConv);           
            
            HITemp11 = IDC.ComputeAreaIntegrationFiniteSupport(shapeGeom,@F11,paramsIJ,doConv);
            HITemp12 = IDC.ComputeAreaIntegrationFiniteSupport(shapeGeom,@F12,paramsIJ,doConv);
            
            HIInt11 = [HITemp11(:,:,1,1), HITemp11(:,:,1,2) ; ...
                       HITemp11(:,:,2,1), HITemp11(:,:,2,2) ];

            HIInt12 = [HITemp12(:,:,1,1), HITemp12(:,:,1,2) ; ...
                       HITemp12(:,:,2,1), HITemp12(:,:,2,2) ];                   
                   
            HIStruct(iS,jS).HIInt11 = HIInt11;
            HIStruct(iS,jS).HIInt12 = HIInt12;
            
            HIStruct(jS,iS).HIInt11 = HIInt11;
            HIStruct(jS,iS).HIInt12 = HIInt12;
        end
    end

    %--------------------------------------------------------------------------

    function z = F11(x,y,params)         
        z = f11(x,y,params);                    
    end

    function z = F12(x,y,params)         
        z = f12(x,y,params);                    
    end

    
%     function z = F11(x,y)         
%         z = f11(x,y,paramsIJ);                    
%     end
% 
%     function z = F12(x,y)         
%         z = f12(x,y,paramsIJ);                    
%     end
    
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