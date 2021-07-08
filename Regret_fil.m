function Filter =  Regret_fil (sys)

sys      = F_NehariLya_in_sys(sys);
N        = F_Nehari_solution(sys);

d        = sys.dimF;
Filter.F = zeros(3*d,3*d);
Filter.H = zeros(d,3*d);
Filter.G = zeros(3*d,d);
Filter.J = zeros(d);

Temp_invhalf    = inv ( chol(sys.RP,'lower'));

F21     = - N.GN*Temp_invhalf*sys.H;
F31     = sys.FW*sys.U*sys.H'*inv(sys.RP)*sys.H ...
         - sys.KQ*N.HN*N.GN*Temp_invhalf*sys.H;
F32     = sys.KQ*N.HN*N.FN;
F1      = [sys.FP,zeros(d),zeros(d)];
F2      = [F21,N.FN,zeros(d)];
F3      = [F31,F32,sys.FW];

Filter.F = [F1;F2;F3];

H1  = sys.L - sys.L*(sys.P-sys.U)*sys.H'*inv(sys.RP)*sys.H...
    - N.HN*N.GN*Temp_invhalf*sys.H;
Filter.H = [H1,N.HN*N.FN,sys.L];


G3  = sys.KQ*N.HN*N.GN*Temp_invhalf...
    -sys.FW*sys.U*sys.H'*inv(sys.RP);
Filter.G = [sys.KP;N.GN*Temp_invhalf;G3];

Filter.J = sys.L*(sys.P-sys.U)*sys.H'*inv(sys.RP) + N.HN*N.GN*Temp_invhalf;


%This checks the minimality
OBS = [Filter.H;Filter.H*Filter.F;Filter.H*Filter.F*Filter.F];

CON = [Filter.G,Filter.F*Filter.G,Filter.F*Filter.F*Filter.G];

svd(OBS);
svd(CON);