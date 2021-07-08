function g_final = Find_gamma(sys)

gmin = 0;
gmax = 5;
run  = 1;

while Check_g(gmax,sys)
gmax = gmax*2;
end

while run
    if gmax-gmin>10e-7
        g_mid = (gmin+gmax)/2;
    else
        g_final = gmax;
        run = 0;
    end
    Fail = Check_g(g_mid,sys);
    if Fail==1 %(if didnt work)
        gmin = g_mid;
    else
        gmax = g_mid;
    end
end



% disp(['gamma = ',num2str(g_final)]);








