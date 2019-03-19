function [Kopen, indop, PsiK] = get_psiK(E, mu, eps, Ls, Q, Ri, Ri1)
//[Kopen, indop, PsiK] = get_psiK(E, mu, eps, Ls, Q, Ri, Ri1)
//
//  Ri  : match point
//  Ri1 : point BEFORE Ri
//

//  GCG, 19-Mar-2003
//       6-Aug-2003 - improve condition for large L channels
//       eps is now channel energy!

Max     = 1d100 // overflow protection
//mprintf('get_psiK: Bessel overflow protection, max = %g\n', Max)

indop   = find(eps <= E)
indcl   = find(eps >  E)
nop     = length(indop)
ncl     = length(indcl)
nbas    = nop+ncl
if Ls==[], Ls=zeros(nbas,1); end

// if nop==1
//   mprintf('%d open channel\n',nop)
// else
//   mprintf('%d open channels\n',nop)
// end

k   = sqrt(abs(2*mu*(E-eps)))
fi  = zeros(nbas,1)
fi1 = zeros(nbas,1)
gi  = zeros(nbas,1)
gi1 = zeros(nbas,1)

// open channeLs
if indop<>[]
  L   = Ls(indop)

  fi(indop)  = s_f(L,k(indop),Ri)
  fi1(indop) = s_f(L,k(indop),Ri1)
  gi(indop)  = s_g(L,k(indop),Ri, Max)  // overflow protection
  gi1(indop) = s_g(L,k(indop),Ri1,Max)  // overflow protection
end

// closed channeLs

if indcl<>[]
  L          = Ls(indcl)
  gi(indcl)  = s_c2(L,k(indcl),Ri, [], Max)     // overflow protection
  gi1(indcl) = s_c2(L,k(indcl),Ri1,Ri, Max)     // overflow protection
end

mu2 = sqrt(mu)
F   = diag(mu2*fi)
F1  = diag(mu2*fi1)
G   = diag(mu2*gi)
G1  = diag(mu2*gi1)
F   = F(:,indop)
F1  = F1(:,indop)

// K   = (G1-Q*G)\(F1-Q*F)
// Use column scaling to improve the condition beyond believe, 6-Aug-2003

B   = G1-Q*G
nrm = zeros(nbas, 1)
for i = 1:nbas
  nrm(i) = norm(B(:,i))
  B(:,i) = B(:,i)/nrm(i)
end
K   = B\(F1-Q*F)
for i = 1:nbas
  K(i, :) = K(i, :)/nrm(i)
end


Kopen = K(indop, :)
if argn(1)>2
  PsiK = F - G*K
end
endfunction

function [f, fp] = s_f(l, k, R)
//S_F   Regular wave
//  [F, FP] = S_F(L, k, R)
//  SQRT(MU)*F is flux normalized regular wave
//  Wronskian: F*GP - GP*F = 1
//
//  See also S_G, S_C1, S_C2

//  GCG, 6-Aug-2003, use SLATEC DBESJ routine
//  GCG, 15-Jul-2010, now use built-in besselj

l    = l(:)
k    = k(:)
R    = R(:)
z    = k.*R
fac  = sqrt(0.5*%pi*R)
//f    = fac.*ff_dbesj(l+0.5, z)
f    = fac.*besselj(l+0.5, z)
if argn(1)>1
//  fp  = (l+1).*f./R - k.*fac.*ff_dbesj(l+1.5, z)
  fp  = (l+1).*f./R - k.*fac.*besselj(l+1.5, z)
end
endfunction

function [g, gp] = s_g(l, k, R, Max)
//S_G  irregular wave
//  [G, GP] = S_G(L, K, R)
//  SQRT(MU)*G is flux normalized irregular wave
//  Wronskian: F*GP - GP*F = 1
//
//  See also S_F, S_C1, S_C2

//  GCG, 6-Aug-2003
//  GCG, 15-Jul-2010, now use built-in bessely
//  GCG, 21-May-2012, protect agains overflow

  nrhs = argn(2)
  if nrhs<4, Max=[], end
  if Max==[], Max=0, end

  l    = l(:)
  k    = k(:)
  R    = R(:)
  z    = k.*R
  fac  = sqrt(0.5*%pi*R)
  [g, ierror]   = evstr('fac.*bessely(l+0.5, z)')
  if ierror<>0
    if Max==0
      error(msprintf('ERROR in s_g: overflow in bessely, ierror = %g', ierror))
    else
      nl   = size(l, '*')
      nz   = size(z, '*')
      nf   = size(fac, '*')
      n    = max([nl nz nf])
      if nl==1
        l    = l*ones(n,1)
      end
      if nz==1
        z    = z*ones(n,1)
      end
      if nf==1
        fac   = fac*ones(n,1)
      end
      g    = Max*ones(n,1)

      for i=1:n
        [g(i), ierr] = evstr('fac(i)*bessely(l(i)+0.5, z(i))')
        if ierr<>0
           mprintf('s_g (lib_cc), WARNING: %d: l=%d: %g*bessely(%g,%g) set to %g\n', i, l(i), fac(i), l(i)+0.5, z(i), g(i))
        end
      end
      ind = find(abs(g)>Max)
      if ind<>[]
        n   = size(ind, '*')
        gt  = sign(g(ind))*Max
        for i=1:n
          mprintf('s_g (lib_cc), WARNING: %d: l=%d: %9.3e truncated to %9.3e\n', ind(i), l(ind(i)),g(ind(i)), gt(i))
        end
        g(ind) = gt
      end
   end
end
if argn(1)>1
//  gp  = (l+1).*g./R - k.*fac.*ff_dbesy(l+1.5, z)
  gp  = (l+1).*g./R - k.*fac.*bessely(l+1.5, z)
end
endfunction

function [c2, c2p] = s_c2(l, k, R, R0, Max)
//S_C2  irregular closed channel (Bessel K)
//  [C2, C2P] = S_C2(L, K, R, R0, Max)
//  SQRT(MU)*C2*EXP(-K*R0) is flux normalized regular closed channel
//  Wronskian: C1*C2P - C1P*C2 = 1
//
//  See also S_C1, S_F, S_G

//  GCG, 6-Aug-2003 
//  GCG, 15-Jul-2010, now use built-in besselk 
//  GCG, 15-Mar-2012, protect besselk against overflow

  [nlhs, nrhs] = argn()

  if nrhs<5, Max  = [], end
  if nrhs<4, R0   = [], end

  if Max==[], Max = 0, end
  if R0==[],  R0  = R($), end

  l            = l(:)
  k            = k(:)
  R            = R(:)
  z            = k.*R
  fac          = sqrt(2*R/%pi)
  scale        = exp(k.*(R0-R))
  [c2, ierror] = evstr('-fac.*besselk(l+0.5, z, 1).*scale')
  if ierror<>0
    if Max==0
      error(msprintf('ERROR in s_g: overflow in besselk, ierror = %g', ierror))
    else

      nl           = size(l, '*')
      nz           = size(z, '*')
      nf           = size(fac, '*')
      n            = max([nl nz nf])
      if nl==1
        l            = l*ones(n,1)
      end
      if nz==1
        z     = z*ones(n,1)
        scale = scale*ones(n,1)
      end
      if nf==1
        fac   = fac*ones(n,1)
      end
      g    = Max*ones(n,1)

      for i=1:n
        [c2(i), ierr] = evstr('-fac(i).*besselk(l(i)+0.5, z(i), 1).*scale(i)')
        if ierr<>0
           mprintf('s_g (lib_cc), WARNING: %d: l=%d: %g*bessely(%g,%g) set to %g\n', i, l(i), fac(i), l(i)+0.5, z(i), g(i))
        end
      end
      ind = find(abs(c2)>Max)
      if ind<>[]
        n   = size(ind, '*')
        c2t = sign(c2(ind))*Max
        for i=1:n
          mprintf('s_c2 (lib_cc), WARNING: %d: l=%d: %9.3e truncated to %9.3e\n', ind(i), l(ind(i)), c2(ind(i)), c2t(i))
        end
        c2(ind) = c2t
      end
    end
  end

  if nlhs>1
    // c2p   = (l+1).*c2./R + k.*fac.*ff_dbesk(l+1.5, z, 2).*scale
    c2p   = (l+1).*c2./R + k.*fac.*besselk(l+1.5, z, 1).*scale
  end
endfunction
