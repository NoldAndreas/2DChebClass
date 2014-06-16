function z = HIweight_11(ptsPolLoc,params)
    ptsCart = Pol2CartPts(ptsPolLoc);
    HIfn = params.HIfn;
    z = HIfn(ptsCart.y1_kv,ptsCart.y2_kv,params);
    z = z(:,1,1);
end