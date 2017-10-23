function vel = velocity_interp_old(x,t,xx, yy,velocityGrid)
    %rentre un vecteur x de taille 2*n avec n le nombre de points
    
    vel = zeros(2,1);
    epsilonx = abs(xx(1,1) - xx(1,2))*0.5;
    epsilony = abs(yy(1,1) - yy(2,1))*0.5;
    
    for i=1:size(xx,1)
        for j=1:size(xx,2)
            deltax = x(1) - xx(i,j);
            deltay = x(2) - yy(i,j);
            if abs(deltax)<epsilonx && abs(deltay)<epsilony
                
                vel = [velocityGrid(i,j,1); velocityGrid(i,j,2)] ;
                
            end
        end
    end

end
