function Pts = shiftPointsCart(pts,y10,y20)
    % only need to shift physical (y) domain           

    Pts.y1_kv = pts.y1_kv + y10;   Pts.y2_kv = pts.y2_kv + y20;

end
