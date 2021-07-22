function [r1, r2, r3, s1, s2, s3, f3] = matrc(k0, dz, iz, rhob, ...
                                  alpw, alpb, ksqw, ksqb, pdu, pdl)
% The tridiagonal matrices.
  np = size(pdu,  2);
  nz = size(rhob, 2);
  nz = nz - 2;
  r1 = zeros(nz + 2, np);
  r2 = zeros(nz + 2, np);
  r3 = zeros(nz + 2, np);
  s1 = zeros(nz + 2, np);
  s2 = zeros(nz + 2, np);
  s3 = zeros(nz + 2, np);
  f1 = zeros(nz + 2,  1);
  f2 = zeros(nz + 2,  1);
  f3 = zeros(nz + 2,  1);

  a1    = k0 ^ 2 / 6;
  a2    = 2 * k0 ^ 2 / 3;
  a3    = k0 ^ 2 / 6;
  cfact = 0.5 / dz ^ 2;
  dfact = 1 / 12;

% New matrices
  id = 1 : iz;
  f1(id) = 1 ./ alpw(id);
  f2(id) = 1;
  f3(id) = alpw(id);
  ksq(id)= ksqw(id);
  id = iz + 1 : nz + 2;
  f1(id) = rhob(id) ./ alpb(id);
  f2(id) = 1 ./ rhob(id);
  f3(id) = alpb(id);
  ksq(id)= ksqb(id);

% Discretization by Galerkin's method.
  i1 = 2; 
  i2 = nz + 1;
  id = i1 : i2;
  c1 =  cfact * f1(id) .* (f2(id-1) +   f2(id))     .* f3(id - 1);
  c2 = -cfact * f1(id) .* (f2(id-1) + 2*f2(id) + f2(id + 1)).* f3(id);
  c3 =  cfact * f1(id) .* (f2(id)   +   f2(id + 1)) .* f3(id + 1);
  d1 =  c1 + dfact * (ksq(id - 1) +     ksq(id)).';
  d2 =  c2 + dfact * (ksq(id - 1) + 6 * ksq(id) + ksq(id + 1)).';
  d3 =  c3 + dfact * (ksq(id)     +     ksq(id + 1)).';

  r1(id, :) = a1 + d1 * pdl;
  r2(id, :) = a2 + d2 * pdl;
  r3(id, :) = a3 + d3 * pdl;
  s1(id, :) = a1 + d1 * pdu;
  s2(id, :) = a2 + d2 * pdu;
  s3(id, :) = a3 + d3 * pdu;

% The matrix decomposition.
  for id = i1 : iz
        rfact = 1 ./ (r2(id, :) - r1(id, :) .* r3(id - 1, :));
        r1(id, :) = r1(id, :) .* rfact;
        r3(id, :) = r3(id, :) .* rfact;
        s1(id, :) = s1(id, :) .* rfact;
        s2(id, :) = s2(id, :) .* rfact;
        s3(id, :) = s3(id, :) .* rfact;
  end

  for id = i2 : -1 : iz+2
        rfact = 1 ./ (r2(id, :) - r3(id, :) .* r1(id + 1, :));
        r1(id, :) = r1(id, :) .* rfact;
        r3(id, :) = r3(id, :) .* rfact;
        s1(id, :) = s1(id, :) .* rfact;
        s2(id, :) = s2(id, :) .* rfact;
        s3(id, :) = s3(id, :) .* rfact;
  end

  r2(iz + 1, :) = r2(iz + 1, :) - r1(iz + 1, :) .* r3(iz, :);
  r2(iz + 1, :) = r2(iz + 1, :) - r3(iz + 1, :) .* r1(iz + 2, :);
  r2(iz + 1, :) = 1 ./ r2(iz + 1, :);

end

