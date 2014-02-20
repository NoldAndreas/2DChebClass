classdef PerformanceTestClass < handle
    properties 
        A;
        b = 2;
        c = 4;
        e = zeros(1200,200);
        Ch
        d = ones(100,2);
        f_pointer
        
    end
        
    methods
        function this = PerformanceTestClass(c)
            this.A = ones(1200,1)*ones(1,1200)*c;            
            this.f_pointer = @this.getCh;
        end
        
        function this = set.Ch(this,v)
            this.Ch = v;           
        end
        
        function c = testCount0(this,N) 
            s = struct('Ch',2*ones(N,1));            
            c = 1;            
            for i = 1:N
                c = c + s.Ch(i);
            end
        end        
        
        function c = testCount(this,N)                         
            c = 1;            
            this.Ch = 2*ones(N,1);
            for i = 1:N
                c = c + this.Ch(i);
            end
        end        
        function c = testCount2(this,N)
            c = 1;
            Chh = 2*ones(N,1);
            for i = 1:N
                c = c + Chh(i);
            end
        end
        
        function x = getCh(this,i)
            x = this.Ch(i);
        end
        
        function x = getChgetCh(this,i)
            x = getCh(this,i);
        end        
        
        function c = sumGetCh(this,N)
            c = 1;
            this.Ch = 2*ones(N,1);
            for i = 1:N
                c = c + getCh(this,i);
            end
        end
        
        function c = sumGetChgetCh(this,N)
            c = 1;
            this.Ch = 2*ones(N,1);
            for i = 1:N
                c = c + getChgetCh(this,i);
            end
        end
        
        function c = sumGetChPointer(this,N)
            c = 1;
            this.Ch = 2*ones(N,1);
            for i = 1:N
                c = c + this.f_pointer(i);
            end
        end
        
    end
    
end