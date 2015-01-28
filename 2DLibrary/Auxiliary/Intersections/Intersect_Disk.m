function dataDisk = Intersect_Disk(MainShape,diskShape)%y20,r,N,sphere)
    r   = diskShape.R;
    N   = [diskShape.N1,diskShape.N2];
    	
	y20 = diskShape.Origin(2);    
    
    if(isa(MainShape,'HalfSpace'))
        %Cartesian shift
        
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end

        shape.N      = N;
        shape.R      = r;     
        shape.Origin = diskShape.Origin;
        if((nargin >= 2) && (diskShape.sphere == true))
            shape.sphere = true;
        end

        %1. find points of disk in HalfSpace            
        if(y20 == inf)
            shape.Origin(2) = 0;
            area            = Disc(shape); 
        elseif(y20 >= y2Min + r)
            %1a. if full disk is in HalfSpace                
            %shape.N         = GetPointNumber(shape,'Disc',this.Accuracy);  
            area            = Disc(shape);                                              
        elseif((y20 < (y2Min + r)) && (y20 >= y2Min))
            %1b. if part of disk is in HalfSpace  (>= half)
            %1b1. Integrate over segment in HalfSpace
            %shape.Origin    = [0,y20];
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Min;       
            
            %shape.h         = y20 - y2Min;
            %shape.Top       = true;

            shape.NW        = [2*N(1),N(2)];
            shape.NT        = N;

            area            = BigSegment(shape);                                                                                
        elseif((y20 < y2Min) && (y20 >= y2Min - r))                
            %shape.h         = y2Min - y20;            
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Min;            
           % shape.N         = GetPointNumber(shape,'Segment',this.Accuracy);
            area            = SegmentWall(shape);                       
        %1c. if part of disk is in HalfSpace  (< half)    
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end            
                
    elseif(isa(MainShape,'InfCapillary'))
        %Cartesian shift
        
        y2Min = MainShape.y2Min;
        y2Max = MainShape.y2Max;
        if(isa(MainShape,'InfCapillarySkewed'))
%            exc = MException('Intersect_Disk:InfCapillarySkewed','case not implemented/tested');
%            throw(exc);                
            y2Min = y2Min*sin(MainShape.alpha);
            y2Max = y2Max*sin(MainShape.alpha);
        end

        shape.N      = N;
        shape.R      = r;     
        shape.Origin = diskShape.Origin;
        if((nargin >= 2) && (diskShape.sphere == true))
            shape.sphere = true;
        end        
        
        if((y20 >= y2Min - r) && (y20 < y2Min))                     
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Min;            
           % shape.N         = GetPointNumber(shape,'Segment',this.Accuracy);
            area            = SegmentWall(shape);      
        elseif(y20 == y2Min)
            shape.Rout = r;
            shape.th1  = 0;
            shape.th2  = pi;
            area = Wedge(shape);                         
        elseif((y20 > y2Min) && (y20 < y2Min + r))
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Min;       
            
            shape.NW        = [2*N(1),N(2)];
            shape.NT        = N;

            area            = BigSegment(shape);
        elseif((y20 >= y2Min + r) && (y20 <= y2Max - r))
            area            = Disc(shape);            
        elseif((y20 > y2Max-r) && (y20 < y2Max))
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Max;       
            
            shape.NW        = [2*N(1),N(2)];
            shape.NT        = N;

            area            = BigSegment(shape);            
        elseif(y20 == y2Max)
            shape.Rout = r;
            shape.th1  = pi;
            shape.th2  = 2*pi;
            area = Wedge(shape);
        elseif((y20 > y2Max) && (y20 <= y2Max+r))                     
            shape.Wall_VertHor = 'horizontal';
            shape.Wall_Y       = y2Max;            
           % shape.N         = GetPointNumber(shape,'Segment',this.Accuracy);
            area            = SegmentWall(shape);                               
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end                    
        
    elseif(isa(MainShape,'Box'))
        area         = Intersect_Disc_Box(diskShape,MainShape);        
    else
        exc = MException('Intersect_Disk','case not implemented');
        throw(exc);                
    end    
    
    dataDisk.pts       = area.GetCartPts();        
    ptsLoc.y1_kv       = dataDisk.pts.y1_kv - area.Origin(1);
    ptsLoc.y2_kv       = dataDisk.pts.y2_kv - area.Origin(2);        
    dataDisk.ptsPolLoc = Cart2PolPts(ptsLoc);
    if(y20 == inf)
        dataDisk.pts.y2_kv = dataDisk.pts.y2_kv + y20;
    end

    [dataDisk.int,dataDisk.area] = area.ComputeIntegrationVector();
        
end   