function convStruct = FexMatrices_ConstShortRange(optsPhys,IDC)
  
    params       = optsPhys.V2;    
        
    shape.N      = optsPhys.FexNum.N;
    shape.volume = true;        
    shape.R      = params.LJsigma;

    conv1 = IDC.ComputeConvolutionFiniteSupport(Sphere(shape),'',IDC.Pts);
    
    shape.R = params.lambda;
    conv2 = IDC.ComputeConvolutionFiniteSupport(Sphere(shape),'',IDC.Pts);

    convStruct.Conv =  -(conv2(:,:,1) - conv1(:,:,1));
    
    %Rescale
    alpha  = -2/3*pi*(params.lambda^3 - params.LJsigma^3);        
    c      = (-16/9*pi)/alpha;    
    convStruct.Conv = c*convStruct.Conv;
    
    if(isfield(params,'epsilon'))
        convStruct.Conv = convStruct.Conv*params.epsilon;
    end

end