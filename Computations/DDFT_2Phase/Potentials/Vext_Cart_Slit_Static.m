function V =Vext_Cart_Slit_Static(y2,optsPhys)
    
    yw_bottom = min(y2); 
    yw_top    = max(y2); 
        
    epsilon_w = optsPhys.epsilon_w; 
    
    d      = yw_top - yw_bottom;
    a2     = AttractiveWall(yw_top-y2,epsilon_w(2),1);
    a4     = AttractiveWall(y2-yw_bottom,epsilon_w(4),1);    
                           
    V      = a2 + a4;                        
                    
end