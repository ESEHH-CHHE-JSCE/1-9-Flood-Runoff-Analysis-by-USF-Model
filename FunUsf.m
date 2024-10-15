%% Urban Storage Function

function F = FunUsf(X,RIEO,k1,k2,k3,p1,p2,z)

s = k1*(X(1)^(p1/p2))+X(2)*k2;

if X(1) > 1.0E-20
    if s >= z 
           F =[X(2) ; -k1/k2*p1/p2*X(2)*X(1)^(p1/p2-1) - 1/k2*X(1)^(1/p2)...
            - k1*k3/k2*X(1)^(p1/p2) - k3*X(2) + 1/k2*(RIEO+ k3*z)];
    else
           F =[X(2) ; -k1/k2*p1/p2*X(2)*X(1)^(p1/p2-1) - 1/k2*X(1)^(1/p2)...
            + 1/k2*(RIEO)];
    end
else
    if s >= z
          F =[X(2) ; -k3*X(2)+(1/k2)*(RIEO+(k3*z))];
    else
          F =[X(2) ; (1/k2)*RIEO];
    end
end 