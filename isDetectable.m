function fail = isDetectable(A,B)

fail=0;
v_eig = eig(A);

for i=1:size(A,1)
    if rank ([v_eig(i)*eye(size(A)) - A; B])<size(A,1)
       fail=1; 
    end
end