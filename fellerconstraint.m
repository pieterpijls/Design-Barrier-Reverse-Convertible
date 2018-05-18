function [C,Ceq] = fellerconstraint(x)
C = x(3)^2-2*x(1)*x(2); % Feller contsraint using Heston parameters to avoid negative vol
Ceq = 0;
end


