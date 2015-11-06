   function y = RotnePragner_f2(pts)
       y = cos(pts.y2_kv).*sin(pts.y2_kv)./(pts.y1_kv).^2;
    end