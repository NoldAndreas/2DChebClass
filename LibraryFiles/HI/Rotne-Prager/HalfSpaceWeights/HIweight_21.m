function z = HIweight_21(ptsPolLoc,params)
    ptsCart = Pol2CartPts(ptsPolLoc);
    HIfn = params.HIfn;
    z = HIfn(ptsCart.y1_kv,ptsCart.y2_kv,params);
    z = z(:,2,1);
end