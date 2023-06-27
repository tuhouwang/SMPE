function [phi, r] = FlatStep(X, Layers, Coll, pade1, pade2, ...
                             np, rs, rt, dr, dep, rho, ...
                             initial, Lowerboundary)
    r  = rs : dr : rt;
    nr = length(r);
    
    % The last two lines of each subblock are the rows of boundary 
    % condition replacement.
    tag = zeros(Layers,1);
    for it = 1 : Layers
       tag(it) = sum(Coll(1:it)+1)-1;
    end

    N = size(X, 1) - 1;
    T = eye (N+1);
    for ip = 1 : np
        A = eye(N+1) + pade1(ip) * X;
        B = eye(N+1) + pade2(ip) * X;
        
        A(tag,  :) = 0.0;
        A(tag+1,:) = 0.0;
        
        % Apply boundary conditions to B.
        B = ImposeCondition(tag, B, dep, rho, Coll, Lowerboundary);       
        T =  B \ A * T;    
    end  
    
    psi = cell(1, nr);
    % phi represents the spectral coefficient stored by layers.
    phi = cell(Layers, nr);
    phi(:,1) = initial;
    psi(  1) = {cell2mat(initial)};
    % Step forward.
    for ir = 2 : nr
        psi(ir) = {T * psi{ir-1}};
        for it = 1 : Layers
           phi{it,ir} = psi{ir}(sum(Coll(1:it-1)+1)+1:sum(Coll(1:it)+1)); 
        end
    end

end
