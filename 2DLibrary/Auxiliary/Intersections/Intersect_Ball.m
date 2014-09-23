function dataBall = Intersect_Ball(MainShape,ballShape)
%function dataBall = Intersect_Ball(MainShape,y20,r,N)
    r   = ballShape.R;
    N   = [ballShape.N1,ballShape.N2];
	y20 = ballShape.Origin(2);       
    
    if(isa(MainShape,'HalfSpace'))   
        
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end
        shape.N  = N;
        shape.R  = r;     

        %1. find points of disk in HalfSpace            
        if(y20 >= y2Min + r)
            %1a. if full disk is in HalfSpace                
            shape.theta1 = 0;
            shape.theta2 = pi;                             
        elseif((y20 < (y2Min + r)) && ...
               (y20 >= (y2Min - r)))
            %1b. if part of disk is in HalfSpace  (>= half)
            %1b1. Integrate over segment in HalfSpace
            %shape.Origin    = [0,y20];
            th                = acos((y20 - y2Min)/r);
            shape.theta1      = 0;
            shape.theta2      = pi-th;                                                
        %1c. if part of disk is in HalfSpace  (< half)    
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end            

        %shape.N            = GetPointNumber(shape,'Ball',MainShape.Accuracy);
        area               = Ball(shape); 
        dataBall.pts       = area.GetCartPts();   %PtsCart;

        %Shift in y2-direction
        dataBall.ptsPolLoc = Cart2PolPts(area.GetCartPts());%PtsCart);            
        dataBall.pts.y2_kv = dataBall.pts.y2_kv + y20;

        [dataBall.int,dataBall.area]     = area.ComputeIntegrationVector();                                   
    elseif(isa(MainShape,'InfCapillary') && (MainShape.y2Max - MainShape.y2Min)> 2*r)   
%         
         y2Min = MainShape.y2Min;
         y2Max = MainShape.y2Max;
         if(isa(MainShape,'InfCapillarySkewed'))   
             y2Min = y2Min*sin(MainShape.alpha);
             y2Max = y2Max*sin(MainShape.alpha);
         end
         shape.N  = N;
         shape.R  = r;     
 
         %1. find points of disk fully in InfCapillary
         if( (y20 >= y2Min + r) && (y20 <= y2Max - r))
%             %1a. if full disk is in HalfSpace                
             shape.theta1 = 0;
             shape.theta2 = pi;                             
         elseif((y20 < (y2Min + r)) && (y20 >= (y2Min - r)))
             %1b. if part of disk is in InfCapillary  (>= half)
             %1b1. Integrate over segment in HalfSpace
             %shape.Origin    = [0,y20];
             th                = acos((y20 - y2Min)/r);
             shape.theta1      = 0;
             shape.theta2      = pi-th;                                                
         elseif((y20 > (y2Max - r)) && (y20 <= (y2Max + r)))
             th                = acos((y20 - y2Max)/r);
             shape.theta1      = th;
             shape.theta2      = pi;  
         %1c. if part of disk is in HalfSpace  (< half)    
         else
             exc = MException('Intersect_Ball','case not implemented');
             throw(exc);                
         end            
% 
         %shape.N            = GetPointNumber(shape,'Ball',MainShape.Accuracy);
         area               = Ball(shape); 
         dataBall.pts       = area.GetCartPts();   %PtsCart;
 
         %Shift in y2-direction
         dataBall.ptsPolLoc = Cart2PolPts(area.GetCartPts());%PtsCart);            
         dataBall.pts.y2_kv = dataBall.pts.y2_kv + y20;
 
         [dataBall.int,dataBall.area]     = area.ComputeIntegrationVector();                                           
    else
        exc = MException('Intersect_Ball','case not implemented');
        throw(exc);                
    end        
end   