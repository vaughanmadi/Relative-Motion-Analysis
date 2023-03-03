%propogating a numerical cartesian state
function[ydot] = yprop(t,y,mu)
rmag = (y(1)^2 + y(2)^2 +y(3)^2)^(1/2);
ydot = zeros(size(y));

ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6);
ydot(4) = -(mu/(rmag^3))*y(1);
ydot(5) = -(mu/(rmag^3))*y(2);
ydot(6) = -(mu/(rmag^3))*y(3);
