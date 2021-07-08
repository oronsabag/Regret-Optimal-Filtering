function Nehari_Const =  F_Nehari_solution (sys)


%This apply to both scenarios
Nehari_Const.GN         = inv( eye(size(sys.FP,1))  - sys.FP*sys.Zg*sys.FP'*sys.PI)*sys.FP*sys.Zg*sys.H'*inv(chol(sys.RP));
Nehari_Const.FN         = sys.FP - Nehari_Const.GN*inv( chol (sys.RP ,'lower'))*sys.H;
Nehari_Const.HN         = sys.L*(sys.P-sys.U)*sys.FP'*sys.PI;

