function psi = KernelFunc(Vec, dz, dep, Layers)

    z   = 0 : dz : dep{end}(end);
    zm  = [];
    psi = [];
    
    for i = 1 : Layers
        zi = dep{i}(1) : dz : dep{i}(end);
        xi = -2 / (dep{i}(end) - dep{i}(1)) * zi + ...
                  (dep{i}(end) + dep{i}(1)) / (dep{i}(end) - dep{i}(1));
        if(isempty(zm) == 1)
           zm = zi;
           psi= InvChebTrans(Vec{i}, xi);
        elseif(zm(end) == zi(1))   
           zm  = [zm, zi(2:end)];
           p   = InvChebTrans(Vec{i}, xi);
           psi = [psi; p(2:end,:)];             
        else
           zm  = [zm, zi];
           psi = [psi; InvChebTrans(Vec{i}, xi)];          
        end  
    end
    
    psi = interp1(zm', psi, z', 'linear', 'extrap');

end