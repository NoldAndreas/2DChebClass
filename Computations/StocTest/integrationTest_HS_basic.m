function integrationTest_HS_basic

R = 1;

params.R = 1;

params.sigmaHS = 0.5;

HSgeom.N = [10;10]; HSgeom.L1 = 2; HSgeom.L2 = 2;
HSgeom.y2wall = 0; HSgeom.N2bound = 10; HSgeom.h = 1;
HSgeom.L2_AD = 1; HSgeom.alpha_deg = 90; HSgeom.R = R;

HS = HalfSpace_FMT(HSgeom);

nPts = length(HS.Pts.y1_kv);

%Ageom.N = [20;20]; Ageom.R = R;
%Ageom.shape = 'Disc';

Ageom.N = [20;20]; Ageom.RMin = R; Ageom.L = 2;
Ageom.shape = 'InfAnnulus';


shapeClass = str2func(Ageom.shape);

Atemp = shapeClass(Ageom);

testFn = str2func('RP12_2D');

if(nargin(testFn)==3)
    useParams = true;
    fPTemp = testFn(Atemp.Pts.y1_kv,Atemp.Pts.y2_kv,params);
else
    useParams = false;
    fPTemp = testFn(Atemp.Pts.y1_kv,Atemp.Pts.y2_kv);
end

fDim = size(fPTemp);

nElts = prod(fDim(2:end));

M_int = zeros([nPts,nPts,fDim(2:end)]);
Mmask = repmat({':'},[1,fDim]);

for iPt = 1:nPts
    Ageom.Origin = [HS.Pts.y1_kv(iPt);HS.Pts.y2_kv(iPt)];
    A = shapeClass(Ageom);
    IntersectArea = Intersect(HS,A);

    Int    = IntersectArea.int;   % 1 x N1*N2  
    IntT   = Int.';               % N1*N2 x 1
    IntT = IntT(:,ones(1,nElts)); % N1*N2 x nElts
    IntT = IntT(:);               % N1*N2*nElts x 1
    
    IP     = HS.SubShapePts(IntersectArea.pts);
    
    if(useParams)
        f = testFn(IntersectArea.pts.y1_kv-HS.Pts.y1_kv(iPt),IntersectArea.pts.y2_kv-HS.Pts.y2_kv(iPt),params);
    else
        f = testFn(IntersectArea.pts.y1_kv-HS.Pts.y1_kv(iPt),IntersectArea.pts.y2_kv-HS.Pts.y2_kv(iPt));
    end
    
    nPtsF = size(f,1);
    
    f = f(:);   % nPtsF*nElts x 1
     
    IntF = (IntT.*f).';   % 1 x nPtsF*nElts
    
    InterpIntF = zeros(nPts,nElts);

    for iElt = 1:nElts
        eltMask = (iElt-1)*nPtsF+1 : iElt*nPtsF;
        InterpIntF(:,iElt) = (IntF(eltMask)*IP).';
    end
     
    InterpIntF = reshape(InterpIntF,[nPts,fDim(2:end)]);
    
    Mmask{1} = iPt;
    M_int(Mmask{:}) = InterpIntF;
    
end


%rho = ones(size(HS.Pts.y1_kv));

rho = exp(-(HS.Pts.y1_kv.^2 + HS.Pts.y2_kv.^2));

val11 = M_int(:,:,1,1)*rho;
val12 = M_int(:,:,1,2)*rho;
val21 = M_int(:,:,2,1)*rho;
val22 = M_int(:,:,2,2)*rho;

[val11 val12 val21 val22 HS.Pts.y1_kv HS.Pts.y2_kv]

    function z = testFa(x,y)
        z = ones(size(x));
    end

    function z = testFb(x,y)
        z = 2*ones(size(x));
    end

    function z = testF2(x,y)
        z = zeros(length(x),2,2);
        z(:,1,1) = testFa(x,y);
        z(:,1,2) = testFb(x,y);
        z(:,2,1) = testFb(x,y);
        z(:,2,2) = testFa(x,y);
    end

end