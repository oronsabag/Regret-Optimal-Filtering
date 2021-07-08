function LyaK = Lyapunov_Hankel(T)

% dlyap(A,B,C); solves AXB-X+C=0;

% F PI F* - pi + G*G = 0

A_PI = T.F;
B_PI = T.F';
C_PI = T.G*T.G';

LyaK.PI = dlyap(A_PI,B_PI,C_PI);

%F'ZF - Z + g^-2H* H = 0;

A_Z = T.F';
B_Z = T.F;
C_Z = T.H'*T.H;

LyaK.Z = dlyap(A_Z,B_Z,C_Z);

