function B = ImposeCondition(tag, B, dep, rho, Coll, Lowerboundary)

        Layers = length(tag);
        %Interface conditions.
        for it = 1 : Layers-1
            B(tag(it), sum(Coll(1:it-1)+1)+1:sum(Coll(1:it)+1)) = (-1.0).^(0:Coll(it));
            B(tag(it), sum(Coll(1:it)+1)+1:sum(Coll(1:it+1)+1)) =  -1.0 ;
            
            D  = DerivationMatrix(Coll(it) + 1);
            Pu = -1 / rho{it}(end) / (dep{it}(end) - dep{it}(1)) * ((-1.0).^(0 : Coll(it))) * D;
            
            D  = DerivationMatrix(Coll(it+1) + 1);
            Pd =  1 / rho{it+1}(1) / (dep{it+1}(end) - dep{it+1}(1)) * ones(1, Coll(it+1)+1) * D;
            B(tag(it)+1, sum(Coll(1:it-1)+1)+1:sum(Coll(1:it)+1)) = Pu;
            B(tag(it)+1, sum(Coll(1:it)+1)+1:sum(Coll(1:it+1)+1)) = Pd;
        end
        %Boundary conditions.
        B(tag(end),         1:Coll(1)+1) = 1.0;
        B(tag(end)+1, end-Coll(end):end) = (-1.0).^(0:Coll(end));
        if(Lowerboundary == 'R')
            D = DerivationMatrix(Coll(end) + 1);
            B(tag(end)+1, end-Coll(end):end) = ...
            B(tag(end)+1, end-Coll(end):end) * D;
        end

end
