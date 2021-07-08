function sys = F_NehariLya_in_sys(sys)

% dlyap(A,B,C); solves AXB-X+C=0;

% F PI F* - pi + G*G = 0

A_PI = sys.FP';
B_PI = sys.FP;
C_PI = sys.H'*inv(eye(size(sys.H,1)) + sys.H*sys.P*sys.H')*sys.H;

sys.PI = dlyap(A_PI,B_PI,C_PI);

%F'ZF - Z + g^-2H* H = 0;

if sys.causal
    
    A_Z = sys.FP;
    B_Z = sys.FP';
    C_Z = sys.FP*(sys.P-sys.U)'*sys.L'*inv(sys.RQ)*sys.L*(sys.P-sys.U)*sys.FP';
    
    sys.Zg = dlyap(A_Z,B_Z,C_Z);
    
else
    
    A_Z = sys.FP;
    B_Z = sys.FP';
    C_Z = (sys.P-sys.U)'*sys.L'*inv(sys.RQ)*sys.L*(sys.P-sys.U);
    
    sys.Zg = dlyap(A_Z,B_Z,C_Z);
    
end