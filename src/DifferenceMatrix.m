function [D,x] = DifferenceMatrix(N)

    x  = cos((0 : N) * pi / N)';
    c  = [2; ones(N-1,1); 2] .* (-1).^(0 : N)';

    X  = repmat(x,1,N+1);
    dX = X - X';

    D = (c * (1 ./ c)') ./ (dX + eye(N + 1));
    D = D - diag(sum(D, 2));
end