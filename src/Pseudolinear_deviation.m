function deviation = Pseudolinear_deviation(Coll, dr, dz, Lowerboundary)
    Create_env(Coll, dr, dz, Lowerboundary );

    clc
% close all
clear 
tic;
% edit 'input_SMPE.txt';

[casename, Layers, np, ns, c0, freq, zs, zr, rmax, dr, depth, ... 
 dz, tlmin, tlmax, N, range, Coll, dep, c, rho, alpha, ...
 Lowerboundary] = ReadEnvParameter('input_SMPE.txt');

if(Lowerboundary == 'A')
    % Acoustic half-space.
    
    % Add another layer of Perfectly Matched Layer with thickness of two 
    % wavelengths under the original layers.
    lambda = 2 * c0 / freq;
    depth  = [depth; depth(end) + lambda];
    Layers = Layers + 1;
    Coll   = [Coll; max(ceil(Coll(end) * 0.25), 20)];
    
    % The sound speed, density and attenuating coefficient in the PML are 
    % consistent with the seabed (the first SSP).
    PMLd = {depth(end-1:end)'};
    PMLc = {[    c{end,1}(end),    c{end,1}(end)]};
    PMLr = {[  rho{end,1}(end),  rho{end,1}(end)]};
    PMLa = {[alpha{end,1}(end),alpha{end,1}(end)]};
    
    PMLd = repmat(PMLd, 1, N);
    PMLc = repmat(PMLc, 1, N);
    PMLr = repmat(PMLr, 1, N);
    PMLa = repmat(PMLa, 1, N);
    
    dep   = [dep;   PMLd];
    c     = [c;     PMLc];
    rho   = [rho;   PMLr];
    alpha = [alpha; PMLa];
end

[w, k0, z, dep, c, rho, alpha] = Initialization(freq, c0, dz, depth, ... 
 Layers, N, Coll, dep, c, rho, alpha);

if(Lowerboundary == 'A')
    % PML uses complex coordinate transformation technique to absorb 
    % downward radiation energy. 
    tau = (dep{end} - depth(end-1)) / lambda;
    tau = 100 * tau .^ 3 ./ (1 + tau .^ 2);
    % The factor g is used to modify the 
    % Depth Operator of the last layer.
    g = 1.0 ./ (1 + tau + 2i * tau );
else
    g = 1;
end
%***************************obtain the initial field***********************
zd = (0 : 0.1 * dz : depth(end))';
[zl, cl] = Column(Layers, dep(:,1), c(:,1));
cw = interp1(zl, cl, zd, 'linear');
[~, ~, ~, ~, ~, starter] = selfstarter(zs, 0.1 * dz, k0, w, cw, ...
                                       np, ns, c0, dr, length(zd));

phi = {};
r   = dr;
for m = 1 : Layers
    % The self-starter is projected onto the spectral space.
    phi(m,1) = {interp1(zd, starter, dep{m,1}, 'linear')};
    phi(m,1) = {ChebTransFFT(Coll(m,1), phi{m,1})};
end

[pade1, pade2] = epade(np, ns, 1, k0, dr);

range = [dr, range(2:end), rmax];
for j = 1 : length(range)-1
    X = DepthOperator(Layers, k0, c0, Coll, dep(:,j), ...
                      c(:,j), rho(:,j), alpha(:,j), g);

    if( range(j+1) - range(j) > dr )
        % Long flat segment.
        [phit, rt] = FlatStep(X, Layers, Coll, pade1, pade2, np, ... 
                              range(j), range(j+1), dr, dep(:,j), ...
                              rho(:,j), phi(:,end), Lowerboundary);

         phi = [phi, phit(:,2:end)];
         r   = [r, rt(2:end)];
    else
         % Long flat section, multiple steps can be calculated without 
         % updating the Depth Operator.     
         phit = OneStep(X, Layers, Coll, pade1, pade2, np, ...
                        dep(:,j), rho(:,j), phi(:,end), Lowerboundary);               

         r   = [r, range(j)];
         phi = [phi, phit];
    end
end

% Inversion of spectral coefficients into physical space.
psi = zeros(length(z), length(r));
for ir = 1 : length(r)
    psi(:,ir) = KernelFunc(phi(:,ir), dz, dep, Layers);
end

psi = psi .* exp(1i * k0 * dr);
u   = psi ./ sqrt(r);

% Cut off the sound field of PML.
if(Lowerboundary == 'A')
    ind = find( z <= depth(end-1) );
    z   = z(ind);
    u   = u(ind,:);  
end


%*********************calculate and  the results***************************
tl_PE    = - 20 * log10( abs(u));

tl_PE = tl_PE(2:100/dz,2:end);
% if(Lowerboundary == 'V')
%     load('analytic_free.mat');
%     if ( z(end)>100 )
%         load('analytic_rigid.mat');
%         tl_PE = tl_PE(2:100/dz,2:end);
%     end
%     if ( z(end)==100 )
%        
%         tl_PE = tl_PE(2:end-1,2:end);
%     end
% end    
% if(Lowerboundary == 'R' )
%     load('analytic_rigid.mat');
%     tl_PE = tl_PE(2:100/dz,2:end);
% end  
tl=Pseudolinear_analyticV(length(r)+1,length(z));

tl=tl(2:end-1,3:end);
% z_a=z_a(2:end-1);
% z=z(2:end-1);
% tl_a = zeros(length(z),length(r));
% for nn = 1:1:length(r)
%     
%     tl_a(:,nn) = interp1(z_a, tl(:,nn), z, 'linear');
% end
%**************************calculate the deviation*************************
deviation = abs(tl-tl_PE);


deviation = sum(sum(deviation))/(length(z)-2) / (length(r)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coll=[10:1:30,40,50,100];
% error=zeros(2,length(Coll))
% error(1,:)=Coll;
% for i=1:1:length(Coll)
%     error(2,i)= Pseudoliner_deviation(Coll(i), 3);
% end



end


