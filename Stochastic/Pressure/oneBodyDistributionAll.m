function oneBodyDistributionAll(stocStruct,optsPhys)

    x = stocStruct.x;
    p = stocStruct.p;
    
    nTimes = size(x,3);
    
    opts.nParticles = optsPhys.nParticles;
    opts.nBinsX = 30;
    opts.nBinsP = 30;
    
    for iTime = 1:10:nTimes
        
        
        xt = x(:,:,iTime);
        pt = p(:,:,iTime);
        
        xt = xt(:);
        pt = pt(:);
        
        oneBodyDistribution(xt,pt,opts);
        
        pause
        
    end
    
    

end