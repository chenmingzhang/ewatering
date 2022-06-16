function n_ele = nod_to_ele(nex,ney,nez,n_nod)
for k=1:nex
    for j  = 1:ney
        for i = 1:nez      
            n_ele(i,j,k) = mean([n_nod(i+1,j,k),n_nod(i+1,j,k+1),n_nod(i+1,j+1,k+1),n_nod(i+1,j+1,k),...
                                    n_nod(i,j,k),n_nod(i,j,k+1),n_nod(i,j+1,k+1),n_nod(i,j+1,k)] )   ;  % use mean method to get the middle of the cell
        end
    end
end
end