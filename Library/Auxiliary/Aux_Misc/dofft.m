function yf = dofft(y,n1,n2)
        yf =  reshape(fft(reshape(y,n2,n1)),n1*n2,1);
end

