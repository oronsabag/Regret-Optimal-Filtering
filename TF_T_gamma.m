function T = TF_T_gamma(sys,gamma)

Ricc_g      = Riccati_gamma(sys,gamma,0);

Const_left  = inv ( chol ( Ricc_g.RQ,'lower'))*sys.L*(sys.P-Ricc_g.U);
Const_right = sys.H'*inv(chol(sys.RP));

if sys.causal
    
    T.H = Const_left*sys.FP';
    
    T.F  = sys.FP';
    
    T.G = Const_right;
    
    T.J = zeros(size(T.H,1),size(T.G,2));
    
else
    
    T.H = Const_left;
    
    T.F  = sys.FP';
    
    T.G = Const_right;
    
    T.J = Const_left*Const_right;
    
end