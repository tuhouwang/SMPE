function [pdu, pdl]=epade( np, ns, ip, k0, dr)
% function [pdu, pdl]=epade( np, ns, ip, k0, c0, dr)
% np number of pade coefficients
% ns stability constraints
% k0 (rad/m)
% c0 (m/s)
% dr range step (m)
%-----------  self starter:ip=0, split-step:ip=1 -----------

  sig=k0*dr;
  n=2*np;

  if ip==1   % split-step Pade approximation
    nu=0;
    alp=0;
  else       % self-starter approximation
    nu=1;
    alp=-0.25;
  end 

% The factorials.
  for ii=1:n
    fact(ii)=factorial(ii);
  end

% The binomial coefficients.
  for ii=1:n+1
    bin(ii,1)=1;
    bin(ii,ii)=1;
  end
  for ii=3:n+1
    for jj=2:ii-1
      bin(ii,jj)=bin(ii-1,jj-1)+bin(ii-1,jj);
    end
  end

% The accuracy constraints.
  [dg,dh1,dh2,dh3]=deriv(n,sig,alp,nu,bin);

  a=zeros(n,n);
  b(1:n)=dg(2:n+1);

  for ii=1:n
    if (2*ii-1<=n) a(ii,2*ii-1)=fact(ii); end
    for jj=1:ii
      if (2*jj<=n) a(ii,2*jj)=-bin(ii+1,jj+1)*fact(jj)*dg(ii-jj+1); end
    end
  end

  jj=1:np;
% The stability constraints.
  if ns>=1
    z1=-3;
    b(n)=-1;
    a(n,2*jj-1)=z1.^jj;
    a(n,2*jj)=0;
  end 

  if ns>=2
    z1=-1.5;
    b(n-1)=-1;
    a(n-1,2*jj-1)=z1.^jj;
    a(n-1,2*jj)=0.0;
  end

  b=b(:);
  b=a\b;

  dh1(1)=1;
  dh1(jj+1)=b(2*jj-1);
  dh2=fndrt(dh1,np);
  pdu(jj)=-1./dh2(jj);

  dh1(1)=1;
  dh1(jj+1)=b(2*jj);
  dh2=fndrt(dh1,np);
  pdl(jj)=-1./dh2(jj);

return

function [dg,dh1,dh2,dh3]=deriv(n,sig,alp,nu,bin)
% The derivatives of the operator function at x=0.

  ci=complex(0,1);

  dh1(1)=0.5*ci*sig;
  dh2(1)=alp;
  dh3(1)=-2.0*nu;
  exp1=-0.5;
  exp2=-1.0;
  exp3=-1.0;

  for ii=2:n
    dh1(ii)=dh1(ii-1)*exp1;
    dh2(ii)=dh2(ii-1)*exp2;
    dh3(ii)=-nu*dh3(ii-1)*exp3;
    exp1=exp1-1.0;
    exp2=exp2-1.0;
    exp3=exp3-1.0;
  end

  dg(1)=1.0;
  dg(2)=dh1(1)+dh2(1)+dh3(1);
  for ii=2:n
    dg(ii+1)=dh1(ii)+dh2(ii)+dh3(ii);
    for jj=1:ii-1
      dg(ii+1)=dg(ii+1)+bin(ii,jj)*(dh1(jj)+dh2(jj)+dh3(jj))*dg(ii-jj+1);
    end
  end

return

function [z]=fndrt(a,n)
% The root-finding subroutine. 

  if n==1
    z(1)=-a(1)/a(2);
    return
  end
  if n~=2

    for k=n:-1:3
      % Obtain an approximate root.
      root=0;
      err=1e-12;
      root=guerre(a,k,root,err,1000);
  
      % Refine the root by iterating five more times.
      err=0;
      root=guerre(a,k,root,err,5);
      z(k)=root;
  
      % Divide out the factor (z-root).
      for ii=k:-1:1
        a(ii)=a(ii)+root*a(ii+1);
      end
      for ii=1:k
        a(ii)=a(ii+1);
      end
    end
  end

% Solve the quadratic equation.
  z(2)=0.5*(-a(2)+sqrt(a(2)^2-4*a(1)*a(3)))/a(3);
  z(1)=0.5*(-a(2)-sqrt(a(2)^2-4*a(1)*a(3)))/a(3);

return

function [z]=guerre(a,n,guess,err,nter)
%     This subroutine finds a root of a polynomial of degree n > 2
%     by Laguerre's method.

  ci=complex(0,1);
  epsb=1e-20;
  z=guess;

% The coefficients of p'(z) and p''(z).

  ii=[1:n]; az=ii.*a(ii+1);
  ii=[1:n-1]; azz=ii.*az(ii+1);

  iter=0; jter=1;
  while 1
    p=a(n)+a(n+1)*z;

    for ii=n-1:-1:1
      p=a(ii)+z*p;
    end
    if abs(p)<epsb return; end

    pz=az(n-1)+az(n)*z;
    for ii=n-2:-1:1
      pz=az(ii)+z*pz;
    end

    pzz=azz(n-2)+azz(n-1)*z;
    for ii=n-3:-1:1
      pzz=azz(ii)+z*pzz;
    end

    % The Laguerre perturbation.
    f=pz/p;
    g=f^2-pzz/p;
    h=sqrt((n-1)*(n*g-f^2));
    amp1=abs(f+h);
    amp2=abs(f-h);
    if amp1>amp2
      dz=-n/(f+h);
    else
      dz=-n/(f-h);
    end

    iter=iter+1;

    % Rotate by 90 degrees to avoid limit cycles. 
    jter=jter+1;
    if jter==10
      jter=1;
      dz=dz*ci;
    end
    z=z+dz;

    if ((abs(dz)>err)&&(iter<nter))
      continue;
    else
      break;
    end
  end

return

