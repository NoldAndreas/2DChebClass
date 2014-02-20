function V = Corner(x,y,epsilon_w)
    Vx = AttractiveWall(x,epsilon_w,1);
    Vy = AttractiveWall(y,epsilon_w,1);
    
    alpha = - epsilon_w*pi^2/2;
    
    V = (Vx + Vy)/2 - (2*alpha - Rectangle(x,y,epsilon_w))/4;
    
end