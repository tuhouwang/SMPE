function [casename, N, np, freq, zs, zr, rmax, dr, H, dz, tlmin, tlmax, ... 
          dep, c, rho, alpha] = ReadEnvParameter(env_file)

    fid           = fopen(env_file);
    casename      = fgetl(fid);
    N             = fscanf(fid, '%d', 1);
    np            = fscanf(fid, '%d', 1);
    freq          = fscanf(fid, '%f', 1);
    zs            = fscanf(fid, '%f', 1);
    zr            = fscanf(fid, '%f', 1);
    rmax          = fscanf(fid, '%f', 1);
    dr            = fscanf(fid, '%f', 1);
    H             = fscanf(fid, '%f', 1);
    dz            = fscanf(fid, '%f', 1);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    n             = fscanf(fid, '%d', 1);
    
    if(H > 0.0 && N > 2 && np > 2)
        Profile  = fscanf(fid, '%f %f', [4, n]);
        dep      = Profile(1, 1:n);
        c        = Profile(2, 1:n);
        rho      = Profile(3, 1:n);
        alpha    = Profile(4, 1:n);
    else
        error('Error! H must greater than 0, N and np must greater than 2 !');
    end
    
    % Check the input underwater sound profile
    if(dep(1) ~= 0.0 || dep(n) ~= H)
        error('Error! input sound profile is unsuitable !');
    end

    if((H / dz - floor(H / dz)) ~= 0 )
        error('Error! The input dz unsuitable !');
    end

    if((rmax / dr - floor(rmax / dr)) ~= 0)
        error('Please reinput the dr and rmax !');
    end

    if(zs <= 0  || zs >= H || zr <= 0 || zr >= H)
        error('zs and zr must be greater than 0 and less than H !');
    end
    
    if(tlmin >= tlmax)
        error('tlmin must less than tlmax !');
    end    

    fclose(fid);

end