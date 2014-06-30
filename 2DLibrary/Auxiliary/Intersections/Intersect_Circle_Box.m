function line = Intersect_Circle_Box(circle,box)
    
    %Initialization
    y10      = circle.Origin(1);
    y20      = circle.Origin(2);
    
    left     = box.y1Min;
    right    = box.y1Max;
    top      = box.y2Max;
    bottom   = box.y2Min;
    
    R = circle.R;

    centralX = ( (y10 >= left + R) && (y10 <= right-R) );
    centralY = ( (y20 >= bottom + R) && (y20 <= top-R) );
    topY     = (y20 <= top + R) && (y20 >= top - R);
    bottomY  = (y20 <= bottom + R) && (y20 >= bottom - R);
    leftX    = (y10 <= left + R) && (y10 >= left - R);
    rightX  = (y10 <= right + R) && (y10 >= right - R);

    topIn  = (y20 <= top) && (y20 > top - R);
    topOut = (y20 <= top + R) && (y20 > top); 

    bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
    bottomOut = (y20 < bottom ) && (y20 >= bottom - R);

    leftIn   = (y10 < left + R) && (y10 >= left);
    leftOut  = (y10 < left) && (y10 >= left - R);

    rightIn  = (y10 <= right) && (y10 > right - R);
    rightOut = (y10 <= right + R) && (y10 > right);

    containsNW = ( (y10-left)^2 + (y20-top)^2 < R^2);
    containsNE = ( (y10-right)^2 + (y20-top)^2 < R^2);
    containsSE = ( (y10-right)^2 + (y20-bottom)^2 < R^2);
    containsSW = ( (y10-left)^2 + (y20-bottom)^2 < R^2);


    Nin  = topIn && centralX;
    Nout = topOut && (centralX || (leftIn && ~containsNW) || (rightIn && ~containsNE) );

    Sin  = bottomIn && centralX;
    Sout = bottomOut && (centralX || (leftIn &&~containsSW) || (rightIn && ~containsSE) );

    Ein  = rightIn && centralY;
    Eout = rightOut && (centralY || (topIn && ~containsNE) || (bottomIn && ~containsSE) );

    Win  = leftIn && centralY;
    Wout = leftOut && (centralY || (topIn && ~containsNW) || (bottomIn && ~containsSW) );


    central = centralX && centralY;

    NWw = containsNW;
    NEw = containsNE;
    SWw = containsSW;
    SEw = containsSE;

    NWwo = leftIn && topIn && ~containsNW;
    NEwo = rightIn && topIn && ~containsNE;
    SWwo = leftIn && bottomIn && ~containsSW;
    SEwo = rightIn && bottomIn && ~containsSE;

    corner = (NWw || NEw || SEw || SWw);

    doubleCut = (NWwo || NEwo || SEwo || SWwo);

    edgeIn = (Nin || Sin || Ein || Win);
    edgeOut = (Nout || Sout || Eout || Wout);

    if(central)
        shape.N = 20;
        shape.R = R;
        shape.Origin = [y10,y20];
        line = Circle(shape);
    elseif(corner)
        shape.N = 20;
        shape.Origin = [y10;y20];
        shape.R = R;

        if(SWw)
            shape.Corner = 'SW';
            shape.CornerPos = [left;bottom];
        elseif(NWw)
            shape.Corner = 'NW';
            shape.CornerPos = [left;top];
        elseif(NEw)
            shape.Corner = 'NE';
            shape.CornerPos = [right;top];
        elseif(SEw)
            shape.Corner = 'SE';
            shape.CornerPos = [right;bottom];
        end

        line = Arc(shape);

    elseif(doubleCut)
        shape.N = 20;
        shape.Origin = [y10;y20];
        shape.R = R;

        if(SWwo)
            shape.Corner = 'SW';
            shape.CornerPos = [left;bottom];
        elseif(NWwo)
            shape.Corner = 'NW';
            shape.CornerPos = [left;top];
        elseif(NEwo)
            shape.Corner = 'NE';
            shape.CornerPos = [right;top];
        elseif(SEwo)
            shape.Corner = 'SE';
            shape.CornerPos = [right;bottom];
        end

        line = DoubleCutCircle(shape);

    elseif(edgeIn || edgeOut)
        shape.N = 20;
        shape.R  = R;
        shape.Origin = [y10;y20];
        if(Sin || Sout)
            shape.h = y20 - bottom;
            shape.WallPos = 'S';
        elseif(Nin || Nout)
            shape.h = top - y20;
            shape.WallPos = 'N';
        elseif(Ein || Eout)
            shape.h = right - y10;
            shape.WallPos = 'E';
        elseif(Win || Wout)
            shape.h = y10 - left;
            shape.WallPos = 'W';
        end

        line = Arc(shape); 

    else

        line = [];
    end
end