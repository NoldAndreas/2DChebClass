function x = InvQuadMapAB(z,L,zMin,zMax)
%sqrt(1/2) -> hS12
%1  -> zMin
%-1 -> zMax
%
% Map1: SqrtMap: [-1,1] -> [-L,L]
% Map2: f      : [-L,L] -> [zMin,zMax],  f(z) = (zMin+zMax)/2 + z
%

    if( ((zMin == -Inf) && (zMax ~= Inf)) || ((zMin ~= -Inf) && (zMax == Inf)))
        err = MException('SqrtMapAB:zMinMaxOutofRange', ...
                    'zmin or zmax ar inf, but not both');
        throw(err);
    end

    LD = (zMax - zMin)/2;    
        
    if( (zMin == -Inf) && (zMax == Inf) )
        z1 = z;
    else
        z1 = z - (zMax + zMin)/2;
    end
    
    x  = InvQuadMap(z1,L/2,LD);

    
     
end