function G = GHR_lambda(t,lambda)

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
    
    G = IP.InterPol*(IntM*OneOverf(y,lambda));       
   
    function f = OneOverf(t,lambda)
        f = (lambda*(t.^2-(sin(t)).^2).*(pi-t+sin(t).*cos(t))+((pi-t).^2-(sin(t)).^2).*(t-sin(t).*cos(t)))./ ...
                (2*sin(t).*(lambda^2*(t.^2-(sin(t)).^2)+2*lambda*(t.*(pi-t)+(sin(t)).^2)+(pi-t).^2-(sin(t)).^2));
        f(t==0)  = 0;
        f(t==pi) = 0;
    end

end