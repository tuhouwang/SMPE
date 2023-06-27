function X = DepthOperator(Layers, k0, c0, Coll, dep, c, rho, alpha, g)

        X  = [];
        for m = 1 : Layers - 1
            n  = (c0 ./ c{m} .* (1.0 + 1i * alpha{m} / (40.0 * pi * ...
                                   log10( exp(1.0) ) ) ) ) .^ 2 - 1.0;

            C  = ConvolutionMatrix( ChebTransFFT(Coll(m), n) );
            D  = DerivationMatrix ( Coll(m) + 1);
            Xm = 4.0 / (dep{m}(end) - dep{m}(1)) ^ 2 / k0 ^ 2 * ...
                 ConvolutionMatrix( ChebTransFFT(Coll(m), rho{m} ) ) * D * ...
                 ConvolutionMatrix( ChebTransFFT(Coll(m), 1.0 ./rho{m}) ) * D + C;
             
            X = blkdiag(X, Xm);           
        end
        % The last layer may be PML, so the Depth Operator shall be 
        % calculated separately.
        m  = Layers;
        n  = (c0 ./ c{m} .* (1.0 + 1i * alpha{m} / (40.0 * pi * ...
                                log10( exp(1.0) ) ) ) ) .^ 2 - 1.0;

        C  = ConvolutionMatrix(ChebTransFFT(Coll(m), n));
        D  = DerivationMatrix (Coll(m) + 1);
        Xm = 4.0 / (dep{m}(end) - dep{m}(1)) ^ 2 / k0 ^ 2 * ...
             ConvolutionMatrix(ChebTransFFT(Coll(m), g.*rho{m})) * D * ...
             ConvolutionMatrix(ChebTransFFT(Coll(m), g./rho{m})) * D + C;
             
        X = blkdiag(X, Xm);       
end
