function vel = velocity_vortex(x,t,positions_vortices, gamma_vortices)
    %rentre un vecteur x de taille 2*n avec n le nombre de points
    
    vel = zeros(2,1);
    delta = 0.1;

    for i=1:length(positions_vortices)
        gamma = gamma_vortices{i};
        temp = positions_vortices{i};
        xc = temp(1);
        yc = temp(2);

        vel(1,:) = vel(1,:) - gamma*(x(2,:)-yc)./(2*pi*((x(1,:)-xc).^2 + (x(2,:)-yc).^2 + delta));
        vel(2,:) = vel(2,:) + gamma*(x(1,:)-xc)./(2*pi*((x(1,:)-xc).^2 + (x(2,:)-yc).^2 + delta));
    end





end

