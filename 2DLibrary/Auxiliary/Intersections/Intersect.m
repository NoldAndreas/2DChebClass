function data = Intersect(MainShape,SecondShape)

    if(isa(SecondShape,'Disc'))
        data = Intersect_Disk(MainShape,SecondShape);
    elseif(isa(SecondShape,'Sphere'))
        data = Intersect_Sphere(MainShape,SecondShape);                
    elseif(isa(SecondShape,'Ball'))
        data = Intersect_Ball(MainShape,SecondShape);
    elseif(isa(SecondShape,'Circle'))
        data = Intersect_Circle(MainShape,SecondShape);
    elseif(isa(SecondShape,'InfAnnulus'))
        data = Intersect_InfAnnulus(MainShape,SecondShape);
    else
        exc = MException('Intersect','case not implemented');
        throw(exc);
    end

end