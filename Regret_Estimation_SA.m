clearvars
close all
clc

sys                         = Patameters();

% Tran_func                   = Transfer_functions(sys);

if sys.causal
    
    Tran_func.RegretF       = Regret_fil(sys);
    
else
    
    Tran_func.RegretF       = Regret_fil_strictly(sys);
    
end