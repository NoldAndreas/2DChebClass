function convStruct = FexMatrices_Meanfield(optsPhys,IDC)
%************************************************************************
%convStruct  = Polar_SpectralFourier_ConvolutionMatrices(Pts,Int,optsPhys)
%DESCRIPTION:
% computes matrix M_conv such that
%           M_conv * g_i = int f(y1-y1t,y2-y2t)*g(y1t,y2t) dy1t dy2t,
%INPUT:
%  f   = function for which convolution is to be computed
%  Pts = data for grid points in computational and physical space:
%        - 'y1_kv' - y1-variable for each point of 2D-grid 
%        - 'y2_kv' - y2-variable for each point of 2D-grid
%        - 'x1'    - grid in computational space only fo 1st variable
%        - 'x2'    - grid in computational space only for 2nd variable
%  Int = Integration weight vector
%OUTPUT:
% M_conv = convolution matrix, size (N1*N2,N1*N2) , see above
%************************************************************************

    if(isfield(optsPhys,'nParticlesS'))
        nSpecies=length(optsPhys.nParticlesS);
    else
        nSpecies=length(optsPhys.nSpecies);
    end

    %paramNames = optsPhys.potParams2Names;
    %nParams    = length(paramNames);
    
    params = optsPhys.V2;
    
    f      = str2func(params.V2DV2);
    
    params  = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams = length(paramNames);
    
%    Phi = @PhiUniversal;
%    if(strcmp(IDC.polar,'polar'))
%         Phi = @PhiPolar;
%    else
%        Phi = @PhiCart;
%    end
    
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
	function z = PhiUniversal(r)         
         z = zeros(length(r),nSpecies,nSpecies);
        
         for iSpecies = 1:nSpecies
             for jSpecies = 1:nSpecies
                 
                optsPhysIJ = getIJParams(iSpecies,jSpecies);

                z(:,iSpecies,jSpecies) = f(r,optsPhysIJ);                
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
