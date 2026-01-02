function phi = OneStep_new(X, Layers, Coll, pade1, pade2, np, ...
                       dep, rho, initial, Lowerboundary)
     
    dP = zeros(1, np);
    for ip = 1 : np
        dP(ip) = (pade1(ip) - pade2(ip)) * prod( (pade1([1:ip-1 ip+1:np])...
                - pade2(ip)) ./ (pade2([1:ip-1  ip+1:np]) - pade2(ip)) );
    end
    
    if(np == 1)
        dP(1) = pade1(1) - pade2(1);
    end
    dP = - dP ./ pade2;
    d0 = 1 - sum(dP);

    % The last two lines of each subblock are the rows of boundary 
    % condition replacement.
    tag = zeros(Layers,1);
    for it = 1 : Layers
       tag(it) = sum(Coll(1:it)+1)-1;
    end
    
    N   = size(X, 1) - 1;
    psi = cell2mat(initial);
    W_total = zeros(N+1, 1);
    for ip = 1 : np
        B = eye(N+1) + pade2(ip) * X;
        A = dP(ip) * psi;
        
        A(tag)   = 0.0;
        A(tag+1) = 0.0;
        
        % Apply boundary conditions to B.
        B = ImposeCondition(tag, B, dep, rho, Coll, Lowerboundary); 
        
        W = B \ A;
        W_total = W_total + W;
    end
    psi = d0 * psi + W_total;
    
    % phi represents the spectral coefficient stored by layers.
    phi = cell(Layers, 1);
    for it = 1 : Layers
        phi{it} = psi(sum(Coll(1:it-1)+1)+1:sum(Coll(1:it)+1)); 
    end

end
