function convStruct = FexMatrices_SplitDisk(optsPhys,IDC)
  
    params = optsPhys.V2;
    
    fstr      = (params.V2DV2);    
    f         = func2str(fstr);
    
    params     = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams    = length(paramNames);
    
    shapeAnn.RMin = params.LJsigma;
    if(isfield(optsPhys,'FexNum'))
        shapeAnn.f    = str2func(fstr);%optsPhys.FexNum.L;
        shapeAnn.N    = optsPhys.FexNum.N;
    else    
        shapeAnn.L    = optsPhys.optsNum.PhysArea.Conv.L;
        shapeAnn.N    = optsPhys.optsNum.PhysArea.Conv.N;
    end
    
    annulusArea       = InfAnnulus(shapeAnn);

    shapeDisc.R = params.LJsigma;
    shapeDisc.N = shapeAnn.N;
    shapeDisc.volume = false;
    %diskArea    = Disc(shapeDisc);
    diskArea    = Sphere(shapeDisc);

    conv1 = IDC.ComputeConvolutionFiniteSupport(annulusArea,{@PhiUniversal},IDC.Pts);
    conv2 = IDC.ComputeConvolutionFiniteSupport(diskArea,{@PhiUniversal},IDC.Pts);

    convStruct.Conv =  conv1(:,:,2) + conv2(:,:,2);          
        
    function z = PhiUniversal(r)         
        z = f(r,params);
    end

end