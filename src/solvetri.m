function upd = solvetri( iz, u, r1, r2, r3, s1, s2, s3)
  
  [nz, np] = size(r1);
  v        = zeros(nz, 1);
  nz       = nz - 2;
  epsb     = 1.0d-30;
  for jj = 1 : np

    % The right side.
    ind = 2 : nz + 1;
    v(ind) = s1(ind, jj) .* u(ind - 1) + s2(ind, jj) .* u(ind) + ...
                                    s3(ind, jj) .* u(ind + 1) + epsb;

    % The elimination steps.
    for ii = 3 : iz
         v(ii) = v(ii) - r1(ii, jj) * v(ii - 1) + epsb;
    end
    for ii = nz : -1 : iz+2
        v(ii) = v(ii) - r3(ii, jj) * v(ii + 1) + epsb;
    end

    u(iz + 1) = (v(iz + 1) - r1(iz + 1, jj) * v(iz) -  r3(iz + 1, jj)...
                                   * v(iz + 2)) * r2(iz + 1, jj) + epsb;

    % The back substitution steps.
    for ii = iz : -1 : 2
        u(ii) = v(ii) - r3(ii, jj) * u(ii + 1) + epsb;
    end
    for ii = iz + 2 : nz + 1
        u(ii) = v(ii) - r1(ii, jj) * u(ii - 1) + epsb;
    end
  end
  upd = u;
end

