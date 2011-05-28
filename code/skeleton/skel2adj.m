function [adj] = skel2adj(skel)

    [r,c,v] = find(skel);
    
    rs = [r(1)];
    cs = [c(1)];
    
    [m n] = size(skel);
    
    adj = sparse(m*n,m*n);
    
    while (length(rs)>0)
        r = rs(1);
        c = cs(1);
        rs = rs(2:end);
        cs = cs(2:end);
        idx = ((r-1)*n)+c;
        if (adj(idx,idx) == 0) 
            adj(idx,idx) = 1;
            if (skel(r+1,c) == 1) 
                idx2 = (r*n)+c;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r+1];
                cs = [cs c];
            end
            if (skel(r-1,c) == 1) 
                idx2 = ((r-2)*n)+c;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r-1];
                cs = [cs c];
            end
            if (skel(r,c+1) == 1) 
                idx2 = (r*n)+c+1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r];
                cs = [cs c+1];
            end
            if (skel(r,c-1) == 1) 
                idx2 = (r*n)+c-1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r];
                cs = [cs c-1];
            end
            if (skel(r+1,c+1) == 1) 
                idx2 = (r*n)+c+1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r+1];
                cs = [cs c+1];
            end
            if (skel(r-1,c+1) == 1) 
                idx2 = ((r-2)*n)+c+1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r-1];
                cs = [cs c+1];
            end
            if (skel(r+1,c-1) == 1) 
                idx2 = (r*n)+c-1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r+1];
                cs = [cs c-1];
            end
            if (skel(r-1,c-1) == 1) 
                idx2 = ((r-2)*n)+c-1;
                adj(idx,idx2) = 1;
                adj(idx2,idx) = 1;
                rs = [rs r-1];
                cs = [cs c-1];
            end
        end
    end
    
end
