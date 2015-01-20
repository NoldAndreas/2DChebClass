function G = GHR_lambdaEta(t,lambdaEta)

    N          = 200;
    shape      = struct('N',N,'yMin',0,'yMax',pi);    
    SL  = SpectralLine(shape);    
 
	IP = SL.ComputeInterpolationMatrix(SL.CompSpace(t));
        
    y    = SL.Pts.y;  
    IntM = zeros(N);
	vh   = zeros(1,N);
    for i = 2:N
        %Matrix for integral int(v(y),y=y..Y)
        hh          = y(i) - y(i-1);
        vh([i-1,i]) = vh([i-1,i]) + hh/2;
        IntM(i,:)    = vh;
    end
    
    G = IP.InterPol*(IntM*OneOverf(y,lambdaEta));       
   
    function f = OneOverf(t,lambdaEta)
        f = 1./f_stokes(t,lambdaEta);
    end

end