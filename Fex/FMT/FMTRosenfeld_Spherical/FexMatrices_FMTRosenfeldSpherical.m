function convStruct = FexMatrices_FMTRosenfeldSpherical(optsPhys,IDC)

    if(isfield(optsPhys,'nParticlesS'))
        nSpecies=length(optsPhys.nParticlesS);
    else
        nSpecies=length(optsPhys.nSpecies);
    end
    
    optsFMT.N = 100;
    
    %preallocate
    convStruct(nSpecies).D0 = [];
    
    sigmaS = diag(optsPhys.V2.sigmaS);
    
    for iSpecies = 1:nSpecies
        %R = optsPhys.sigmaS(iSpecies)/2;
        R = sigmaS(iSpecies)/2;
        optsFMT.yMin = -R;
        optsFMT.yMax = R;
        % f2, f3, f2_r, f2r_w
        convTemp = IDC.ComputeConvolutionMatrixPointwise(@collectedF,optsFMT);
        Mf2    = convTemp(:,:,1);
        Mf3    = convTemp(:,:,2);
        Mf2_r  = convTemp(:,:,3);
        Mf2r_w = convTemp(:,:,4);
        convStruct(iSpecies).D2   = Mf2;
        convStruct(iSpecies).D3   = Mf3;
        convStruct(iSpecies).D2r  = Mf2_r;
        convStruct(iSpecies).D0   = 1/(4*pi*R^2)*Mf2;
        convStruct(iSpecies).D1   = 1/(4*pi*R)*Mf2;
        convStruct(iSpecies).D1r  = 1/(4*pi*R)*Mf2_r;
        convStruct(iSpecies).D2rF = Mf2r_w;
        convStruct(iSpecies).D1rF = 1/(4*pi*R)*Mf2r_w;
    end
    

       
%--------------------------------------------------------------------------


    %----------------------------------------------------------------------
    % Functions
    %----------------------------------------------------------------------

    % note that r should be a scalar and s a vector
    
    function f = collectedF(r,s)
        f(:,:,1) = f2(r,s);
        f(:,:,2) = f3(r,s);
        f(:,:,3) = f2_r(r,s);
        f(:,:,4) = f2r_w(r,s);
    end
    
    function [f] = f2(r,s)
        % integral is 1/r \int ds s w2(r-s) \rho(s)
        f = bsxfun(@times,1./r,s) .* weight_w2( bsxfun(@plus, r ,  -s ));
    end

    function [f] = f3(r,s) 
        % integral is 1/r \int ds s w3(r-s) \rho(s)
        f =  bsxfun(@times,1./ r,s) .* weight_w3( bsxfun(@plus, r ,  -s ));
    end

    function [f] = f2_r(r,s) 
        % integral is 1/r \int ds s W2z(r-s) \rho(s) + 1/r^2 \int ds s W3(r-s) \rho(s)
        f =  bsxfun(@times,1./ r,s) .* weight_W2z( bsxfun(@plus, r ,  -s ))+ ...
            + bsxfun(@times,1./(r.^2),s) .* weight_w3( bsxfun(@plus, r ,  -s ));
    end

    function [f] = f2r_w(r,s) 
        % this is for the integral \int \delta N_2(s)/\delta \rho(r)
        % \partial \Phi / \partial N_2 where N_2(r) = n_3(r)/r e_r 
        % + 1/r \int ds \rho(s) s W_2(r-s).
        f  =  bsxfun(@times,1./r,s) .* weight_W2z( bsxfun(@plus,-r ,  s ))+ ...
            + bsxfun(@times,1./(r.^2),s) .* weight_w3( bsxfun(@plus, -r ,  s ));
    end

    %**********************************************
    % Weight functions
    %**********************************************

    function w = weight_w2(r)
        % weight function is 2 \pi R \Theta(R-r)
        r = abs(r);
        markl1 = (r < R);
        w = 2*markl1*pi*R;
        
    end
    
    function w = weight_w3(r)
        % weight function is \pi (R^2 - r^2) \Theta(R-r)
        r = abs(r);
        markl1 = (r < R);
        w =  pi*(R^2 - r.^2).*markl1;  
        
    end

    function w = weight_W2z(r)
        % weight function is 2 \pi r \Theta(R-r)
        markl1 = (abs(r) < R);
        w =  2*pi*r.*markl1;
    end


end
