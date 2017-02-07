function yf = doifft(y,n1,n2)
        yf =  reshape(ifft(reshape(y,n2,n1)),n1*n2,1);
end

