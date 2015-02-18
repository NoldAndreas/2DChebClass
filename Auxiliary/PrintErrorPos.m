function [error,ind,cnt] = PrintErrorPos(error,text,Pts1,Pts2)    

    cnt = 0;

    if(sum(isnan(error))>0)
        cnt = cnt + fprintf(['Error of ',text,': ']);
        if(sum(isnan(error))==1)
            cnt = cnt + cprintf('*red','There is a NaN at ');
            ind = isnan(error);
            if(nargin==3)
                cnt = cnt + fprintf([' (y) = (',num2str(Pts1(ind)),')\n']);    
            elseif(nargin==4)
                cnt = cnt + fprintf(['  (y1,y2) = (',num2str(Pts1(ind)),...
                                 ' , ',num2str(Pts2(ind)),')\n']);    
            end
        else
            cnt = cnt + cprintf('*red','There is a NaN in more that one position in the vector!\n');
        end
        return;
    end

    [error,ind] = max(abs(error));
    cnt = cnt + fprintf(['Error of ',text,': ']);
    if(error < 1e-3)
        str2 = num2str(error,'%10.2e');             
    else
        str2 = num2str(error);     
    end
  
     if(error > 10^-3)
         cnt = cnt + cprintf('*red',str2);
     elseif(error <= 10^(-8))
         cnt = cnt + cprintf('*blue',str2);
     else
         cnt = cnt + cprintf('*magenta',str2);
     end
     pause(0.01) %this is because otherwise, the colors are mixed (I dont know why, but this seems to do the job).
     if(nargin == 3)         
         if(isstruct(Pts1))
            cnt = cnt + fprintf([' , max at (y1,y2) = (',num2str(Pts1.y1_kv(ind)),...
                                 ' , ',num2str(Pts1.y2_kv(ind)),')\n']);    
         else
             cnt = cnt + fprintf([' , max at y = ',num2str(Pts1(ind)),'\n']);    
         end
     elseif(nargin==4)
         cnt = cnt + fprintf([' , max at (y1,y2) = (',num2str(Pts1(ind)),...
                                 ' , ',num2str(Pts2(ind)),')\n']);    
     else
         cnt = cnt + fprintf('.');
         cnt = cnt + fprintf('\n');
     end

end