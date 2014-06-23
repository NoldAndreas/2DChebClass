function params = JeffreyOnishiPreprocess(params)

    % these should be nSpecies x nSpecies
    sigma  = params.sigmaHS;
    
    sigmaII = diag(sigma);
    sigmaII = sigmaII(:);
    
    nSpecies = length(sigmaII);
    
    sigma1  = repmat(sigmaII,1,nSpecies);
    sigma2  = sigma1.';
    
    lambda  = max(sigma1,sigma2)./min(sigma1,sigma2);
        
    params.lambda = lambda;
    
    params.yMin  = sigma; 
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
    
	params.fx = fx;
    
    params.nMax = 11;
   
end