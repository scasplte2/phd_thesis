function W = wmat_sj(R,pars)

W=2*pars.mu*(exp(-pars.alpha*R)*pars.Vymat+diag(pars.lambda));

end
