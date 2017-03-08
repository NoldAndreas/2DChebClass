function x = InvTrefInterval(y,y1,y2,y0,ep_y)

    d  = InvLinearMap(y0,y1,y2);
    ep = 2*ep_y/(y2-y1);
    
    x2       = InvLinearMap(y,y1,y2);
    x        = InvTref(x2,d,ep);
    

end