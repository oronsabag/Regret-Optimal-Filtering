function valid = sys_dim(sys)

valid = 1;
if size(sys.F,1)*size(sys.L,2)*size(sys.H,2)*size(sys.G,1) - size(sys.F,2)^4
    
    valid = 0;
    
end