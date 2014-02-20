function x = InvLinearMap(z,a,b)        
        % Finite linear mapping from [a,b] -> [-1,1]
        x = -1 + 2*(z-a)/(b-a);        
end