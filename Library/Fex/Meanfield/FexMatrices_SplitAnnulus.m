function convStruct = FexMatrices_SplitAnnulus(optsPhys,IDC)        

    params     = optsPhys.V2;
        
    f          = str2func(params.V2DV2);
    
    params     = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams    = length(paramNames);
    
    shapeAnn.RMin = params.LJsigma;
    shapeAnn.RMax = params.r_cutoff;
    if(isfield(optsPhys,'FexNum'))        
        shapeAnn.N    = optsPhys.FexNum.N;
    else            
        shapeAnn.N    = optsPhys.optsNum.PhysArea.Conv.N;
    end
       
    shapeDisc.R = params.LJsigma;
    shapeDisc.N = shapeAnn.N;
    shapeDisc.volume = false;
    %diskArea    = Disc(shapeDisc);

    annulusArea = Annulus(shapeAnn);
    %conv1       = IDC.ComputeConvolutionFiniteSupport(annulusArea,{fstr},IDC.Pts);
    conv1       = IDC.ComputeConvolutionFiniteSupport(annulusArea,{@PhiUniversal},IDC.Pts);
    diskArea    = Sphere(shapeDisc);
    conv2       = IDC.ComputeConvolutionFiniteSupport(diskArea,{@PhiUniversal},IDC.Pts);
    
    %shapeDisc.R = r_cutoff;
    %diskArea    = Sphere(shapeDisc);
    %conv3      = IDC.ComputeConvolutionFiniteSupport(diskArea,{fstr},IDC.Pts);

    convStruct.Conv =  conv1(:,:,2) + conv2(:,:,2); %conv3(:,:,2);
        
    function z = PhiUniversal(r)         
        z = f(r,params);
    end


end