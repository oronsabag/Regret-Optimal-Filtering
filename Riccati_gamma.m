function output = Riccati_gamma(sys,gamma,I)

%   [X,L,G] = DARE(A,B,Q,R,S,E) solves
%                                      -1 
%    E'XE = A'XA - (A'XB + S)(B'XB + R)  (A'XB + S)' + Q
%          G = (B'XB + R)  (B'XA + S'),

% First factorization
% W = sys.F' W sys.F+ QW  - KW' RW KW 

QW = sys.H'*sys.H + gamma^(-2)*sys.L'*sys.L;

[output.W,eig_W,output.KW] = dare(sys.F,sys.G,QW,eye(size(sys.G,2)));
output.FW = sys.F - sys.G*output.KW;

RW =  eye(size(sys.G,2)) + sys.G'*output.W*sys.G;


% [output.Q,eig_Temp,Temp_KQ] = dare(output.FW',Lg',Temp,eye(size(output.Lg,1)));
% 
% output.KQ = Temp_KQ';
% output.FQ = output.FW - output.KQ*Lg;
% output.RQ = gamma^2*eye(size(sys.L,1)) + sys.L*output.Q*sys.L';

[output.Q,eig_Temp,Temp_KQ] = dare(output.FW',sys.L',-sys.G*inv(RW)*sys.G',gamma^2*eye(size(sys.L,1)));

output.KQ = Temp_KQ';
output.FQ = output.FW - output.KQ*sys.L;
output.RQ = gamma^2*eye(size(sys.L,1)) + sys.L*output.Q*sys.L';

output.U = dlyap(output.FQ, sys.FP', output.KQ*sys.L*sys.P*sys.FP');

if I
    sys.W = output.W;
    sys.KW = output.KW;
    sys.FW = output.FW;
    sys.RW = RW;
    sys.Q = output.Q;
    sys.KQ = output.KQ;
    sys.FQ = output.FQ;
    sys.U = output.U;
    sys.RQ = output.RQ;
    output = sys;
end