function phi = OneStep(X, Layers, Coll, pade1, pade2, np, ...
                       dep, rho, initial, Lowerboundary)
                     
    % The last two lines of each subblock are the rows of boundary 
    % condition replacement.
    tag = zeros(Layers,1);
    for it = 1 : Layers
       tag(it) = sum(Coll(1:it)+1)-1;
    end
    
    N   = size(X, 1) - 1;
    psi = cell2mat(initial);
    for ip = 1 : np
        A = eye(N+1) + pade1(ip) * X;
        B = eye(N+1) + pade2(ip) * X;
        
        psi = A * psi;
        
        psi(tag)   = 0.0;
        psi(tag+1) = 0.0;
        
        % Apply boundary conditions to B.
        B = ImposeCondition(tag, B, dep, rho, Coll, Lowerboundary); 
        
        psi =  B \ psi;
    end  
    
    % phi represents the spectral coefficient stored by layers.
    phi = cell(Layers, 1);
    for it = 1 : Layers
        phi{it} = psi(sum(Coll(1:it-1)+1)+1:sum(Coll(1:it)+1)); 
    end

end
