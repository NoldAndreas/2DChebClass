function y = HI_Weight_Test(pts)

    r            = pts.y1_kv;
    
    weight       = 1 - exp(-(r./(1-r)).^2);
    weight(r>1)  = 1;
    y            = log(r).*weight;
    y(r == inf)  = 0;
    y(r == 0)    = 0;

end