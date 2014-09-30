function [error,ind] = PrintErrorPos(error,text,Pts1,Pts2)    
    
    if(sum(isnan(error))>0)
        fprintf(['Error of ',text,': ']);
        if(sum(isnan(error))==1)
            cprintf('*red','There is a NaN at ');
            ind = isnan(error);
            if(nargin==3)
                fprintf([' (y) = (',num2str(Pts1(ind)),')\n']);    
            elseif(nargin==4)
                fprintf(['  (y1,y2) = (',num2str(Pts1(ind)),...
                                 ' , ',num2str(Pts2(ind)),')\n']);    
            end
        else
            cprintf('*red','There is a NaN in more that one position in the vector!\n');
        end
        return;
    end

    [error,ind] = max(abs(error));
    fprintf(['Error of ',text,': ']);
    str2 = num2str(error,'%10.2e');     
  
     if(error > 10^-1)
         cprintf('*red',str2);
     elseif(error <= 10^(-5))
         cprintf('*blue',str2);
     else
         cprintf('*magenta',str2);
     end
     pause(0.01) %this is because otherwise, the colors are mixed (I dont know why, but this seems to do the job).
     if(nargin == 3)         
         if(isstruct(Pts1))
            fprintf([' , max at (y1,y2) = (',num2str(Pts1.y1_kv(ind)),...
                                 ' , ',num2str(Pts1.y2_kv(ind)),')\n']);    
         else
             fprintf([' , max at y = ',num2str(Pts1(ind)),'\n']);    
         end
     elseif(nargin==4)
         fprintf([' , max at (y1,y2) = (',num2str(Pts1(ind)),...
                                 ' , ',num2str(Pts2(ind)),')\n']);    
     else
         fprintf('.');
         fprintf('\n');
     end

end