function [dataCircle] = Intersect_Circle(MainShape,circle)            

    y20     = circle.Origin(2);
    r       = circle.R;
    Ncircle = circle.N;

    if(isa(MainShape,'HalfSpace'))
        
        y2Min = MainShape.y2Min;
        if(isa(MainShape,'HalfSpaceSkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
        end
        
        shapeLine.N = Ncircle;
        shapeLine.R = r;            

        %Get Shape for Line
        if(y20 >= y2Min + r)
            line                 = Circle(shapeLine);                                                
        elseif((y20 < (y2Min + r)) && ...
               (y20 >= (y2Min - r)))
            th                 = acos((y20 - y2Min)/r);
            shapeLine.th1      = -pi/2 + th;
            shapeLine.th2      = 3/2*pi-th;

            line               = Arc(shapeLine);
            %dataCircle.pts     = Pol2CartPts(line.Pts);
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end

        dataCircle.pts       = Pol2CartPts(line.Pts);   
        dataCircle.pts.y2_kv = dataCircle.pts.y2_kv + y20;
        dataCircle.ptsPolLoc = line.Pts;            

        [dataCircle.int,dataCircle.length] = line.ComputeIntegrationVector();                                                                         
    elseif(isa(MainShape,'Box'))
        line           = Intersect_Circle_Box(circle,MainShape);        
        dataCircle.pts = line.GetCartPts();
        
        [dataCircle.ptsPolLoc.y2_kv,...
         dataCircle.ptsPolLoc.y1_kv] = cart2pol(dataCircle.pts.y1_kv-circle.Origin(1),...
                                              dataCircle.pts.y2_kv-circle.Origin(2));
     
        [dataCircle.int,dataCircle.length] = line.ComputeIntegrationVector();        
        
        
    elseif(isa(MainShape,'InfCapillary') && (MainShape.y2Max - MainShape.y2Min)>= 2*r)   
%         
        y2Min = MainShape.y2Min;
        y2Max = MainShape.y2Max;
        if(isa(MainShape,'InfCapillarySkewed'))   
            y2Min = y2Min*sin(MainShape.alpha);
            y2Max = y2Max*sin(MainShape.alpha);
        end   
 
        shapeLine.N = Ncircle;
        shapeLine.R = r;
        
        if((y20 < (y2Min + r)) && (y20 >= (y2Min - r)))
            shapeLine.h        = y20 - y2Min;
            shapeLine.WallPos  = 'S';
            line               = Arc(shapeLine);
        elseif( (y20 >= y2Min + r) && (y20 <= y2Max - r))           
            line               = Circle(shapeLine);
        elseif((y20 > (y2Max - r)) && (y20 <= (y2Max + r)))
            shapeLine.h        = y2Max - y20;
            shapeLine.WallPos  = 'N';
            line               = Arc(shapeLine);
        end

        dataCircle.pts       = Pol2CartPts(line.Pts);   
        dataCircle.pts.y2_kv = dataCircle.pts.y2_kv + y20;
        dataCircle.ptsPolLoc = line.Pts;            

        [dataCircle.int,dataCircle.length] = line.ComputeIntegrationVector(); 
         
    else
        exc = MException('HalfSpace_FMT:Intersect_Circle','case not implemented');
        throw(exc);                
    end
end    