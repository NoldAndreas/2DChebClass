function z = GHR_lambda0(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        z = real(z);
end