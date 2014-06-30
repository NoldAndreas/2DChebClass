function data = Intersect(MainShape,SecondShape,opts)

    if(nargin < 3)
        opts = [];
    end

    if(isa(SecondShape,'Disc'))
        data = Intersect_Disk(MainShape,SecondShape,opts);
    elseif(isa(SecondShape,'Ball'))
        data = Intersect_Ball(MainShape,SecondShape,opts);
    elseif(isa(SecondShape,'Circle'))
        data = Intersect_Circle(MainShape,SecondShape,opts);
    elseif(isa(SecondShape,'InfAnnulus'))
        data = Intersect_InfAnnulus(MainShape,SecondShape,opts);
    else
        exc = MException('Intersect','case not implemented');
        throw(exc);
    end

end