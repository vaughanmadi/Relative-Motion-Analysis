% this function is for converting from orbit elements to cartesian 
function [rvec,vvec] = orbtocart(a,e,i,W,w,t,mu,sing)%sing means singularity variable u or l

rPQW = (a*(1-(e^2)))/(1+e*cos(t));
vPQW = sqrt(mu/(a*(1-(e^2))));

if e~=0 && i~=0 %regular degular
rvecPQW = [rPQW*cos(t) rPQW*sin(t) 0];
vvecPQW = [-vPQW*sin(t) vPQW*(e+cos(t)) 0];
Qijk = [cos(W)*cos(w)-cos(i)*sin(W)*sin(w) cos(w)*sin(W)+cos(W)*cos(i)*sin(w) sin(w)*sin(i);-cos(W)*sin(w)-cos(w)*cos(i)*sin(W) cos(W)*cos(w)*cos(i)-sin(W)*sin(w) cos(w)*sin(i);sin(W)*sin(i) -cos(W)*sin(i) cos(i)];
Q = Qijk';
end

if e~=0 && i==0
rvecPQW = [rPQW*cos(t) rPQW*sin(t) 0];
vvecPQW = [-vPQW*sin(t) vPQW*(e+cos(t)) 0];
Qijk = [cos(sing) sin(sing) 0; -sin(sing) cos(sing) 0; 0 0 1]; %sing is wbar
Q = Qijk';
end

if (e==0 && i~=0)
   rvecPQW = [rPQW*cos(sing) rPQW*sin(sing) 0]; % sing = u
   vvecPQW = [-vPQW*sin(sing) vPQW*(e+cos(sing)) 0];
   Qijk = [cos(W) sin(W) 0; -cos(i)*sin(W) cos(W)*cos(i) sin(i);sin(W)*sin(i) -cos(W)*sin(i) cos(i)];
   Q = Qijk';
end

if (e==0 && i==0)
    rvecPQW = [rPQW*cos(sing) rPQW*sin(sing) 0]; %sing = l
    vvecPQW = [-vPQW*sin(sing) vPQW*(e+cos(sing)) 0];
    Q = [1 0 0; 0 1 0; 0 0 1];
end
 rvec = Q*rvecPQW';
 vvec = Q*vvecPQW';
