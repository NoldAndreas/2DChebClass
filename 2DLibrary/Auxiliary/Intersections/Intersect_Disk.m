function dataDisk = Intersect_Disk(MainShape,diskShape,opts)%y20,r,N,sphere)
    r   = diskShape.R;
    N   = [diskShape.N1,diskShape.N2];
    
    
    if(isa(MainShape,'HalfSpace'))
        %Cartesian shift
        y20 = opts.offset_y2;

        
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end                
        
        shape.N  = N;
        shape.R  = r;     
        if((nargin >= 2) && (diskShape.sphere == true))
            shape.sphere = true;
        end

        %1. find points of disk in HalfSpace            
        if(y20 >= y2Min + r)
            %1a. if full disk is in HalfSpace                
            %shape.N         = GetPointNumber(shape,'Disc',this.Accuracy);  
            area            = Disc(shape);                                              
        elseif((y20 < (y2Min + r)) && (y20 >= y2Min))
            %1b. if part of disk is in HalfSpace  (>= half)
            %1b1. Integrate over segment in HalfSpace
            %shape.Origin    = [0,y20];
            shape.h         = y20 - y2Min;
            shape.Top       = true;

            shape.NW        = [2*N(1),N(2)];
            shape.NT        = N;

            area            = BigSegment(shape);                                                                                
        elseif((y20 < y2Min) && (y20 >= y2Min - r))                
            shape.h         = y2Min - y20;                
           % shape.N         = GetPointNumber(shape,'Segment',this.Accuracy);
            area            = Segment(shape);                       
        %1c. if part of disk is in HalfSpace  (< half)    
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end            

        %Shift in y2-direction
        dataDisk.pts       = area.GetCartPts();
        %dataDisk.ptsPolLoc = Cart2PolPts(area.Pts);            
        dataDisk.ptsPolLoc = Cart2PolPts(dataDisk.pts);
        dataDisk.pts.y2_kv = dataDisk.pts.y2_kv + y20;

        [dataDisk.int,dataDisk.area] = area.ComputeIntegrationVector();
        
    elseif(isa(MainShape,'Box'))
        
    else
        exc = MException('Intersect_Disk','case not implemented');
        throw(exc);                
    end
end   