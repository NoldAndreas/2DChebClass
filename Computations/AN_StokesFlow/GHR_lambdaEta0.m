function z = GHR_lambdaEta0(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        if(max(imag(z)) > 1e-12)
            disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        end
        z = real(z);
end