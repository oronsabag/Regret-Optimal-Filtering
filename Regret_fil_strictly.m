function Filter =  Regret_fil_strictly (sys)

sys      = F_NehariLya_in_sys(sys);
N        = F_Nehari_solution(sys);

d        = sys.dimF;
Filter.F = zeros(3*d,3*d);
Filter.H = zeros(d,3*d);
Filter.G = zeros(3*d,d);

RP_invhalf    = inv ( chol(sys.RP,'lower'));

F21     = - N.GN*RP_invhalf*sys.H;
F31     = (sys.FQ*sys.U + sys.KQ*sys.L*sys.P)*sys.H'*inv(sys.RP)*sys.H;
F32     = sys.KQ*sys.L*(sys.P-sys.U)*sys.PI;
F1      = [sys.FP,zeros(d),zeros(d)];
F2      = [F21,N.FN,zeros(d)];
F3      = [F31,F32,sys.FW];

Filter.F = [F1;F2;F3];

Filter.H = [sys.L,sys.L*(sys.P-sys.U)*sys.PI,sys.L];


G3  = - ( sys.FQ*sys.U + sys.KQ*sys.L*sys.P)*sys.H'*inv(sys.RP);
Filter.G = [sys.KP;N.GN*RP_invhalf;G3];

Filter.J = zeros(size(Filter.H,1),size(Filter.G,2));


%This checks the minimality
OBS = [Filter.H;Filter.H*Filter.F;Filter.H*Filter.F*Filter.F];

CON = [Filter.G,Filter.F*Filter.G,Filter.F*Filter.F*Filter.G];

svd(OBS);
svd(CON);