function V =Vext_Cart_Capillary_Static(y1,y2,optsPhys)
    
    yw_left   = min(y1); 
    yw_right  = max(y1); 
    yw_bottom = min(y2); 
    yw_top    = max(y2); 
        
    epsilon_w = optsPhys.epsilon_w; 
    
    d      = yw_top - yw_bottom;
    a1     = Slit(y2-yw_bottom,yw_right-y1,d,epsilon_w(1));
    a2     = AttractiveWall(yw_top-y2,epsilon_w(2),1);
    a3     = Slit(y2-yw_bottom,y1 - yw_left,d,epsilon_w(3));    
    a4     = AttractiveWall(y2-yw_bottom,epsilon_w(4),1);    
                           
    V      = a1 + a2 + a3 + a4;                        
                    
end