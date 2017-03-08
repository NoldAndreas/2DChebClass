function x = mmod(m,n)
   x    = mod(m,n);
   if(x == 0)
       x = n;
   end
end