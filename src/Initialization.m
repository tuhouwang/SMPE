function [w, k0, z, deps, cs, rhos, alphas] = Initialization(freq, ... 
          c0, dz, depth, Layers, N, Coll, dep, c, rho, alpha)

    w  = 2 * pi * freq;
    k0 = w / c0;
    z  = 0 : dz : depth(end);

    deps   = cell(Layers, N);
    cs     = cell(Layers, N);
    rhos   = cell(Layers, N);
    alphas = cell(Layers, N);
    for j = 1 : N
        for m = 1 : Layers
            x           = cos( (0 : Coll(m)) * pi / Coll(m) )';  
            deps{m,j}   = - 0.5 * (dep{m,j}(end) - dep{m,j}(1)) * x ... 
                          + 0.5 * (dep{m,j}(end) + dep{m,j}(1)); 

            cs{m,j}     = interp1(dep{m,j},     c{m,j}, deps{m,j}, 'linear');
            rhos{m,j}   = interp1(dep{m,j},   rho{m,j}, deps{m,j}, 'linear');
            alphas{m,j} = interp1(dep{m,j}, alpha{m,j}, deps{m,j}, 'linear');     
        end
    end
end