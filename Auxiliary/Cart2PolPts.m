function ptsPol = Cart2PolPts(ptsCart)
    ptsPol  = ptsCart;
    [ptsPol.y2_kv,ptsPol.y1_kv] = cart2pol(ptsCart.y1_kv,ptsCart.y2_kv);
end