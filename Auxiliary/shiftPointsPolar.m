function Pts = shiftPointsPolar(pts,y10,y20,polar)
    % only need to shift physical (y) domain           
    %pts are given in polar coordinates!!
    %the fourth variable 'polar' refers to y10,y20

    if(strcmp(polar,'polar'))
          Y10 = y10*cos(y20);    Y20 = y10*sin(y20);
          Y20((y20 == 0) && (y10 == inf)) = 0;
          Y10((y20 == 0) && (y10 == inf)) = inf;
          
          Y20((y20 == pi) && (y10 == inf)) = 0;
          Y10((y20 == pi) && (y10 == inf)) = -inf;
          
          Y20((y20 == pi/2) && (y10 == inf)) = inf;
          Y10((y20 == pi/2) && (y10 == inf)) = 0;
          
          Y20((y20 == 3*pi/2) && (y10 == inf)) = -inf;
          Y10((y20 == 3*pi/2) && (y10 == inf)) = 0;          
          
          %Y10(isnan(Y10)) = 0;   Y20(isnan(Y20)) = 0;        
    else
        Y10 = y10;  Y20 = y20;
    end    

    r = pts.y1_kv; t = pts.y2_kv;

    X = Y10 + r.*cos(t); Y = Y20 + r.*sin(t);

    if(strcmp(polar,'polar'))
        [Pts.y2_kv,Pts.y1_kv] = cart2pol(X,Y); 
        
        mask            = ((r==inf) | (r==-inf));
        Pts.y2_kv(mask) = t(mask);
        Pts.y1_kv(mask) = r(mask);
        
    else
        Pts.y1_kv = X;   Pts.y2_kv = Y;
    end
end
