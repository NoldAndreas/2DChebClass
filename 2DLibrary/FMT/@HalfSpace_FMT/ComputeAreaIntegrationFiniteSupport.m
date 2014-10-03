function M_int = ComputeAreaIntegrationFiniteSupport(this,areaGeom,f,params,convolution)
    shapeClass = str2func(areaGeom.shape);

    nPts = length(this.Pts.y1_kv);

    areaTemp = shapeClass(areaGeom);
    
    doConv  = (nargin > 4 && convolution);
    doParams = ( nargin > 3 && ~isempty(params) );
    
    if(doParams || ~doConv)
        params.pts.x0 = 0;
        params.pts.y0 = 0;
        fPTemp = f(areaTemp.Pts.y1_kv,areaTemp.Pts.y2_kv,params);
    else
        fPTemp = f(areaTemp.Pts.y1_kv,areaTemp.Pts.y2_kv);
    end

    fDim = size(fPTemp);

    nElts = prod(fDim(2:end));

    M_int = zeros([nPts,nPts,fDim(2:end)]);
    Mmask = repmat({':'},[1,fDim]);

    hw = waitbar(0,'Computing Integration Matrix...');
    
    for iPt = 1:nPts
        waitbar(iPt/nPts,hw);
        
        areaGeom.Origin = [this.Pts.y1_kv(iPt) ; this.Pts.y2_kv(iPt)];
        area       = shapeClass(areaGeom);

        IntersectArea = Intersect(this,area);

        Int    = IntersectArea.int;     % 1 x N1*N2  
        IntT   = Int.';                 % N1*N2 x 1
        IntT   = IntT(:,ones(1,nElts)); % N1*N2 x nElts
        IntT   = IntT(:);               % N1*N2*nElts x 1

        IP     = this.SubShapePts(IntersectArea.pts);
        pts    = IntersectArea.pts;
        
        if(doConv)
            pts.y1_kv = pts.y1_kv - this.Pts.y1_kv(iPt);
            pts.y2_kv = pts.y2_kv - this.Pts.y2_kv(iPt);
        else
            params.pts.x0 = this.Pts.y1_kv(iPt);
            params.pts.y0 = this.Pts.y2_kv(iPt);
        end
         
        if(doParams || ~doConv)
            fP = f(pts.y1_kv,pts.y2_kv,params);
        else
            fP = f(pts.y1_kv,pts.y2_kv);
        end
  
        nPtsF = size(fP,1);

        fP = fP(:);   % nPtsF*nElts x 1

        IntF = (IntT.*fP).';   % 1 x nPtsF*nElts

        InterpIntF = zeros(nPts,nElts);

        for iElt = 1:nElts
            eltMask = (iElt-1)*nPtsF+1 : iElt*nPtsF;
            InterpIntF(:,iElt) = (IntF(eltMask)*IP).';
        end

        InterpIntF = reshape(InterpIntF,[nPts,fDim(2:end)]);
        
        InterpIntF(~isfinite(InterpIntF)) = 0;

        Mmask{1} = iPt;
        M_int(Mmask{:}) = InterpIntF;

    end

    delete(hw);

end