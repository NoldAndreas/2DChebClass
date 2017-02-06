function x = InvM1Tref(x2,d,ep)

    c1  = (asinh((1-d)./ep)+asinh((1+d)./ep))/2;
    c2  = asinh((1-d)./ep);

    x   = (asinh((x2-d)./ep)-c2)./c1+1;
end