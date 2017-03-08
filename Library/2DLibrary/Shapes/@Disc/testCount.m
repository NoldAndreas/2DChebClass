function c = testCount(this,N)                         
    c = 1;
    this.Ch = 2*ones(N,1);
    for i = 1:N
        c = c + this.Ch(i);
    end
end        