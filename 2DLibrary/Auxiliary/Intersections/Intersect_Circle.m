function [dataCircle] = Intersect_Circle(MainShape,y20,r,Ncircle)            

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
            dataCircle.pts     = Pol2CartPts(line.Pts);
        else
            exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
            throw(exc);                
        end

        dataCircle.pts       = Pol2CartPts(line.Pts);   
        dataCircle.pts.y2_kv = dataCircle.pts.y2_kv + y20;
        dataCircle.ptsPolLoc = line.Pts;            

        [dataCircle.int,dataCircle.area] = line.ComputeIntegrationVector();                                                                         
        
	else
        exc = MException('Intersect_Circle','case not implemented');
        throw(exc);                
    end
end    