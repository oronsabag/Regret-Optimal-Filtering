function  Fail = Check_g(g,sys)


    T = TF_T_gamma(sys,g);
    Fail = Check_Hankel(T);
