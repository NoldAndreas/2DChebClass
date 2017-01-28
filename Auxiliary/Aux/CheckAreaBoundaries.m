function CheckAreaBoundaries(PhysArea,Area)

    areaName = inputname(2);
    str1 = ['CheckAreaBoundaries:',areaName,'out of range'];

    %Area Limits
    if((Area.y1Max < PhysArea.y1Min) || (Area.y1Max > PhysArea.y1Max))
         err = MException(str1,'y1Max value is outside expected range');
         throw(err);
    end

    if((Area.y1Min < PhysArea.y1Min) || (Area.y1Min > PhysArea.y1Max))
         err = MException(str1,'y1Min value is outside expected range');
         throw(err);
    end
    if((Area.y2Max < PhysArea.y2Min) || (Area.y2Max > PhysArea.y2Max))
         err = MException(str1,'y2Max value is outside expected range');
         throw(err);
    end
    
    if((Area.y2Min < PhysArea.y2Min) || (Area.y2Min > PhysArea.y2Max))
         err = MException(str1,'y2Min value is outside expected range');
         throw(err);
    end

end