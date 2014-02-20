function invertedPts = invertPts(Pts,polar)

    invertedPts = Pts;

    if(strcmp(polar,'polar'))
        tInv = mod(Pts.y2_kv + pi,2*pi);
        invertedPts.y2_kv = tInv;
    else
        xInv = -Pts.y1_kv;
        yInv = -Pts.y2_kv;
        invertedPts.y1_kv = xInv;
        invertedPts.y2_kv = yInv;
    end


end

