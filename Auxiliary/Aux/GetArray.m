    function z = GetArray(min,max,n)
        z = min + (max-min)*(0:(n-1))'/(n-1);
    end