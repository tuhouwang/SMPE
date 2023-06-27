function [casename, Layers, np, ns, c0, freq, zs, zr, rmax, dr, depth, ...
          dz, tlmin, tlmax, N, range, Coll, dep, c, rho, alpha, ...
          Lowerboundary] = ReadEnvParameter(env_file)

    fid           = fopen(env_file);
    casename      = fgetl(fid);
    Layers        = fscanf(fid, '%d', 1);
    np            = fscanf(fid, '%d', 1);
    ns            = fscanf(fid, '%d', 1);
    c0            = fscanf(fid, '%f', 1);
    freq          = fscanf(fid, '%f', 1);
    zs            = fscanf(fid, '%f', 1);
    zr            = fscanf(fid, '%f', 1);
    rmax          = fscanf(fid, '%f', 1);
    dr            = fscanf(fid, '%f', 1);
    depth         = fscanf(fid, '%f', Layers);
    dz            = fscanf(fid, '%f', 1);
    Coll          = fscanf(fid, '%f', Layers);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    N             = fscanf(fid, '%d', 1);
    depth      =  [0; depth];
    range      =  zeros(1, N);
    dep        =  cell(Layers, N);
    c          =  cell(Layers, N);  
    rho        =  cell(Layers, N);
    alpha      =  cell(Layers, N);
    n          =  zeros(Layers);
    
    for j = 1 : N
          range(j) = fscanf(fid, '%f', 1);
          for m = 1 : Layers
              n(m) = fscanf(fid, '%d', 1);
          end
          for m = 1 : Layers
              Profile      = fscanf(fid, '%f %f', [4, n(m)]);
              dep  (m, j)  = {Profile(1, 1:n(m))};
              c    (m, j)  = {Profile(2, 1:n(m))};
              rho  (m, j)  = {Profile(3, 1:n(m))};
              alpha(m, j)  = {Profile(4, 1:n(m))};
              
              %≈–∂œ≤‚…Ó «∑Ò∆•≈‰£°
              if( dep{m,j}(1) ~= depth(m) || dep{m,j}(end) ~= depth(m+1))
                    error('Error! bathymetry and interfaces do not match!');
              end          
          end

    end
    
    Lowerboundary = fscanf(fid, '%s', 1);
    if (Lowerboundary ~= 'V' && Lowerboundary ~= 'R' && Lowerboundary ~= 'A')
        error('Error! The lower boundary must be vaccum, rigid or halfspace!');
    end
    
    if( Layers < 1)
        error('At least one layers of media!');
    end
    
    if( np < 2)
        error('np must greater than or equal to 2!');
    end
    
    depth = depth(2:end);
    if(depth <= 0.0)
        error('Error! H must greater than 0!');
    end
    
    if(N < 1)
        error('At least one bathymetrics are required to judge the terrain!');
    end
    
    if(rmax <= range(end))
        error('Error! rmax must greater than last range!');
    end
    
    if(min(Coll) < 5)
        error('The truncation order is too small to guarantee accuracy!');
    end
    
    % Check the input underwater sound profile
    if(zs <= 0  || zs >= depth(end) || zr <= 0 || zr >= depth(end))
        error('The source and receiver must be in the media!');
    end
    
    if(tlmin >= tlmax)
        error('tlmin must less than tlmax!');
    end    

    fclose(fid);

end