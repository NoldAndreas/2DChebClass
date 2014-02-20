function ptsCart = Pol2CartPts(ptsPol)

    r = ptsPol.y1_kv;
    t = ptsPol.y2_kv;

    [y1,y2] = pol2cart(t,r);

    % fix problematic points
    y2((t == 0) & (r == inf)) = 0;
    y1((t == 0) & (r == inf)) = inf;

    y2((t == pi) & (r == inf)) = 0;
    y1((t == pi) & (r == inf)) = -inf;

    y2((t == pi/2) & (r == inf)) = inf;
    y1((t == pi/2) & (r == inf)) = 0;

    y2((t == 3*pi/2) & (r == inf)) = -inf;
    y1((t == 3*pi/2) & (r == inf)) = 0;

    ptsCart = ptsPol;
    ptsCart.y1_kv = y1;
    ptsCart.y2_kv = y2;
    
end