% The outputs is the system and its dimensions

function sys = Patameters()

sys.causal = 0 ; % 1 stands for causal and 0 is for strictly causal


sys.F = 0.9 * eye(1);
% sys.F  = diag(vec(rand(1,2)));
        sys.F  = rand(2);
% %   sys.F  = rand_stable(3); %This function generates random stable matrix 
% % % sys.F = 
% T = 0.01;
% sys.F = 0.9*eye(2);
% sys.F(1,2) = T;
% sys.G = [0;T];
% sys.H = [1,0];
% sys.L=sys.H;

sys.L       = eye(size(sys.F));
sys.H       = eye(size(sys.F));
sys.G       = eye(size(sys.F));
% 
sys.H       = rand(randi(4,1),size(sys.F,1));
sys.L       = rand(randi(4,1),size(sys.F,1));
sys.G       = rand(size(sys.F,1),randi(4,1));
%  
% sys.H       = rand(randi(4,1),size(sys.F,1));
% sys.L       = rand(1,size(sys.F,1));
% sys.G       = rand(size(sys.F,1),3);
% sys.L=sys.H;
sys.dimF = size(sys.F,1);

if ~sys_dim(sys)
    disp('Dimensions are not compatible');
end

if isDetectable(sys.F,sys.H)
    disp('System not detectable.');
end


[sys.P,eig_P,K_Temp] = dare(sys.F',sys.H',sys.G*sys.G',eye(size(sys.H,1)));
sys.KP = K_Temp';
sys.FP = sys.F - sys.KP*sys.H;
sys.RP = eye(size(sys.H,1)) + sys.H*sys.P*sys.H';

sys.g = Find_gamma(sys);

disp(['Optimal gamma = ',num2str(sys.g)]);

disp(['Regret (squared gamma) = ',num2str(sys.g^2)]);

sys =  Riccati_gamma(sys,sys.g,1); %The '1' is for the first time
