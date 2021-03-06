function HIStruct = HIMatrices2D_AD(opts,IDC)

    optsPhys = opts.optsPhys;
    optsNum  = opts.optsNum;

    if(isfield(optsPhys,'nParticlesS'))
        nSpecies=length(optsPhys.nParticlesS);
    else
        nSpecies=optsPhys.nSpecies;
    end

    params = optsPhys.HI;
    optsNum = optsNum.HINum;
    
    isfield(params,'HIPreprocess')
    
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
    
    for iS = 1:nSpecies
        for jS = iS:nSpecies
            paramsIJ = getIJParams(iS,jS);
            
            paramsIJ.R = paramsIJ.sigmaHS;
    
            HITemp11_Disc = IDC.ComputeConvolutionMatrix(@F11,paramsIJ);  
            HITemp12_Disc = IDC.ComputeConvolutionMatrix(@F12,paramsIJ); 
          
            paramsIJ = rmfield(paramsIJ,'R'); 
            paramsIJ.L = 2;
            paramsIJ.RMin = paramsIJ.sigmaHS;
            paramsIJ.RMax = 50;
            
            HITemp11_Full = IDC.ComputeConvolutionMatrix(@F11,paramsIJ);  
            HITemp12_Full = IDC.ComputeConvolutionMatrix(@F12,paramsIJ); 

%             paramsIJ.RMin = paramsIJ.sigmaHS;          
%             HITemp11_Annulus = IDC.ComputeConvolutionMatrix(@F11,paramsIJ); 
%             HITemp12_Annulus = IDC.ComputeConvolutionMatrix(@F12,paramsIJ); 
            
            HITemp11 = HITemp11_Full - HITemp11_Disc;
            HITemp12 = HITemp12_Full - HITemp12_Disc;
                        
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
    function z = F11(x,y)         
        paramsIJ.HIfn = f11;
        %z = f11(x,y,paramsIJ);
        z = MollifyHIFunction(x,y,paramsIJ);
    end

    function z = F12(x,y) 
        paramsIJ.HIfn = f12;
        %z = f12(x,y,paramsIJ);
        z = MollifyHIFunction(x,y,paramsIJ);
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