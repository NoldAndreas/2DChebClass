function z = ApplyMollifier(x,y1,y2,moll)
    %x: dimension 1: samples 
    %   dimension 2: x1,y1,x2,y2,x3,y3,... where where i is the particle        
    
    z = zeros(length(y1),1);    
    
    for i = 1:length(y1)
        for j = 1:(size(x,2)/2)
            z(i) = z(i) + sum(moll(sqrt( (x(:,2*j-1) - y1(i)).^2 + (x(:,2*j) - y2(i)).^2 )));
        end
    end
    
    z = z/size(x,1);
    
   
    
    
end