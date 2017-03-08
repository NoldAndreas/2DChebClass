function params = JeffreyOnishiPreprocessSpherical(params)

    % these should be nSpecies x nSpecies
    sigma  = params.sigmaHS;
    
    sigmaII = diag(sigma);
    sigmaII = sigmaII(:);
    
    nSpecies = length(sigmaII);
    
    sigma1  = repmat(sigmaII,1,nSpecies);
    sigma2  = sigma1.';
    
    lambda  = max(sigma1,sigma2)./min(sigma1,sigma2);
        
    params.lambda = lambda;
    
    params.yMin  = 0; 
    params.yMax = 50;
    
    fx=zeros(nSpecies,nSpecies,11);
    fx(:,:,1)=3*lambda;
    fx(:,:,2)=9*lambda;
    fx(:,:,3)=-4*lambda + 27*lambda.^2 - 4*lambda.^3;
    fx(:,:,4)=-24*lambda + 81*lambda.^2 + 36*lambda.^3;
    fx(:,:,5)=72*lambda.^2 + 243*lambda.^3 + 72*lambda.^4;
    fx(:,:,6)=16*lambda + 108*lambda.^2 + 281*lambda.^3 + 648*lambda.^4 + 144*lambda.^5;
    fx(:,:,7)=288*lambda.^2 + 1620*lambda.^3 + 1515*lambda.^4 + 1620*lambda.^5 + 288*lambda.^6;
    fx(:,:,8)=576*lambda.^2 + 4848*lambda.^3 + 5409*lambda.^4 + 4524*lambda.^5 ...
       + 3888*lambda.^6 + 576*lambda.^7;
    fx(:,:,9)=1152*lambda.^2 + 9072*lambda.^3+ 14752*lambda.^4+ 26163*lambda.^5 ...
       + 14752*lambda.^6 + 9072*lambda.^7 + 1152*lambda.^8;
    fx(:,:,10)=2304*lambda.^2 + 20736*lambda.^3 + 42804*lambda.^4 + 115849*lambda.^5 ...
       + 76176*lambda.^6 + 39264*lambda.^7 + 20736*lambda.^8 + 2304*lambda.^9;
    fx(:,:,11)=4608*lambda.^2 + 46656*lambda.^3 + 108912*lambda.^4 + 269100*lambda.^5 ...
       + 319899*lambda.^6 + 269100*lambda.^7 + 108912*lambda.^8 ...
       + 46656*lambda.^9 + 4608*lambda.^10;
    
    fy=zeros(nSpecies,nSpecies,11);
    fy(:,:,1)=3/2*lambda;
    fy(:,:,2)=9/4*lambda;
    fy(:,:,3)=2*lambda + 27/8*lambda^2 + 2*lambda^3;
    fy(:,:,4)=6*lambda + 81/16*lambda^2 + 18*lambda^3;
    fy(:,:,5)=63/2*lambda^2 + 243/32*lambda^3 + 63/2*lambda^4;
    fy(:,:,6)=4*lambda + 54*lambda^2 + 1241/64*lambda^3 + 81*lambda^4 + 72*lambda^5;
    fy(:,:,7)=144*lambda^2 + 1053/8*lambda^3 + 19083/128*lambda^4 + 1053/8*lambda^5 ...
           + 144*lambda^6;
    fy(:,:,8)=279*lambda^2 + 4261/8*lambda^3 + 126369/256*lambda^4 - 117/8*lambda^5 ...
           + 648*lambda^6 + 288*lambda^7;
    fy(:,:,9)=576*lambda^2 + 1134*lambda^3 + 60443/32*lambda^4 + 766179/512*lambda^5 ...
           + 60443/32*lambda^6 + 1134*lambda^7 + 576*lambda^8;
    fy(:,:,10)=1152*lambda^2 + 7857/4*lambda^3 + 98487/16*lambda^4 + 10548393/1024*lambda^5 ...
           + 67617/8*lambda^6 - 351/2*lambda^7 + 3888*lambda^8 + 1152*lambda^9;
    fy(:,:,11)=2304*lambda^2 + 7128*lambda^3 + 22071/2*lambda^4 + 2744505/128*lambda^5 ...
           + 95203835/2048*lambda^6 + 2744505/128*lambda^7 + 22071/2*lambda^8 ...
           + 7128*lambda^9 + 2304*lambda^10;
   
	params.fx = fx;
    params.fy = fy;
    
    params.nMax = 11;
   
end