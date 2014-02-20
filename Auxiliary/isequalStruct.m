    function eq = isequalStruct(s1,s2)
        
        eq       = true;        
        f_names1 = fieldnames(s1);
        tol      = 1e-12;
        
        for j = 1:length(f_names1)
            fld = f_names1{j};
            if(~isfield(s2,fld))
              eq = false; return;
            else
                if(isnumeric(s1.(fld)))
                    if(abs(s1.(fld)- s2.(fld)) > tol)
                        eq = false; return;
                    end
                else
                    if(~strcmp(s1.(fld),s2.(fld)))
                        eq = false; return;
                    end
                end
            end            
        end
       
    end