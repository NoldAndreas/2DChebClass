function dataIntersect = Intersect_InfAnnulus(MainShape,infAnnulusShape)
    %Input: 
    %opts: struct with offset_y2
    
    r   = infAnnulusShape.RMin;
    N   = [infAnnulusShape.N1,infAnnulusShape.N2];
    
    y20 = infAnnulusShape.Origin(2);
    
    if(isa(MainShape,'HalfSpace'))
        
        %Cartesian shift     
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end                
        
        shape.N  = N;
        shape.R  = r;     

        %1. find points of disk in HalfSpace            
        if(y20 < y2Min)
            error('Intersect_InfAnnulus: Case not yet implemented');  
        elseif(y20 == Inf)        
            area            = infAnnulusShape;            
            
            dataIntersect.pts       = area.GetCartPts();
            dataIntersect.pts.y2_kv(:) = inf;
            dataIntersect.ptsPolLoc = area.Pts;            
        else
            shape.Origin    = infAnnulusShape.Origin;
            shape.Origin(2) = y20;
            shape.L1        = infAnnulusShape.L;
            shape.y2Wall    = y2Min;            
            area            = HalfSpaceMinusDisk(shape);
                        
            dataIntersect.pts       = area.GetCartPts();
            dataIntersect.ptsPolLoc = area.GetPts();
        end            

        %Shift in y2-direction

        [dataIntersect.int] = area.ComputeIntegrationVector();
        dataIntersect.area  = Inf;
        dataIntersect.shape = area;
    elseif(isa(MainShape,'InfCapillary'))
        
        %Cartesian shift     
        y2Min = MainShape.y2Min;
        y2Max = MainShape.y2Max;
        if(isa(MainShape,'InfCapillarySkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end                
        
        shape.N  = N;
        shape.R  = r;     

        %1. find points of disk in HalfSpace            
        if((y20 < y2Min) || (y20 > y2Max))
            error('Intersect_InfCapillary: Case not yet implemented');          
        else
            shape.Origin    = infAnnulusShape.Origin;
            shape.Origin(2) = y20;
            shape.L1        = infAnnulusShape.L;
            shape.y2Min     = y2Min;            
            shape.y2Max     = y2Max;
            area            = InfCapillaryMinusDisk(shape);
                        
            dataIntersect.pts       = area.GetCartPts();
            dataIntersect.ptsPolLoc = area.GetPts();
        end            

        %Shift in y2-direction
        dataIntersect.shape = area;
        [dataIntersect.int] = area.ComputeIntegrationVector();
        dataIntersect.area  = Inf;
    else
        exc = MException('Intersect_InfAnnulus ','case not implemented');
        throw(exc);                
    end
end   