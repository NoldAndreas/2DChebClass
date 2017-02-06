function [y1t,y2t] = InvSkewedGrid(y1,y2,alpha)
    if(alpha==pi/2)
        y1t = y1;
    else        
        y1t = (y1 - y2/tan(alpha));
    end
    y2t = y2/sin(alpha);
    
    y1t(y1==inf)  = inf;
    y1t(y1==-inf) = -inf;
end
