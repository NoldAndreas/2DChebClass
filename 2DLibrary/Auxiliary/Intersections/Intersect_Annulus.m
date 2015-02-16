function dataIntersect = Intersect_Annulus(MainShape,AnnulusShape)    
    y20 = AnnulusShape.Origin(2);    
    
    if(isa(MainShape,'HalfSpace'))
        
        %Cartesian shift     
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end                        
        h     = y20 - y2Min;
        
        shape.N      = [AnnulusShape.N1,AnnulusShape.N2];
        shape.RMax   = AnnulusShape.RMax;
        shape.RMin   = AnnulusShape.RMin;
        shape.Origin = AnnulusShape.Origin;        

        if(y20 < y2Min)            
            error('Intersect_Annulus: Case not yet implemented');  
        elseif(h >  shape.RMax)        
            area            = AnnulusShape;                                    
        else
            shape.h         = h;
            shape.R_in      = shape.RMin;
            shape.R_out     = shape.RMax;
            area            = AnnulusCut(shape);                                    
        end            
        
        dataIntersect.pts  = area.GetCartPts();            
                
        ptsLoc.y1_kv       = dataIntersect.pts.y1_kv - shape.Origin(1);
        ptsLoc.y2_kv       = dataIntersect.pts.y2_kv - shape.Origin(2);        
        dataIntersect.ptsPolLoc = Cart2PolPts(ptsLoc);       
        
        
        [dataIntersect.int,dataIntersect.area] = area.ComputeIntegrationVector();          
        dataIntersect.shape     = area;    
    else
        exc = MException('Intersect_InfAnnulus ','case not implemented');
        throw(exc);                
    end
end   