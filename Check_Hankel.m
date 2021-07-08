function Fail = Check_Hankel(T)

LyaK = Lyapunov_Hankel(T);

Fail = (max(svd(LyaK.Z*LyaK.PI))>=1);
