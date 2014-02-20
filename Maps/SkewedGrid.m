function [y1,y2] = SkewedGrid(y1t,y2t,alpha)   
    if(alpha == pi/2)
        y1 = y1t;
    else
        y1 = y1t + y2t*cos(alpha);
    end
    y2 = y2t*sin(alpha);
    y1(y1t==inf)  = inf;
    y1(y1t==-inf) = -inf;
end
