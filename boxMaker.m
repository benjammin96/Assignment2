function [cmap] = boxMaker(cmap,box,sigma,nx,ny)

    for j=1:nx
        for i=1:ny     
                if(abs(i-box(1,2)) < box(1,3)/2 & abs(j-box(1,1)) < box(1,4)/2)  
                    cmap(i,j)=sigma;   
                end
        end
    end
end
