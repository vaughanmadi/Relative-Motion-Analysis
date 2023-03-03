clc
clear
%% Environment Analysis
mu = 398600;
elements = [8000 0.0 deg2rad(85) 0 0 0];
p_RSW = [0.025 0 0]'; %in km
p_dot_RSW = [0.001 0 0]'; % in km/s
[r,v] = orbtocart(elements(1),elements(2),elements(3),0,0,0,mu,0);
a = elements(1);
R = r/norm(r);
W = cross(r,v)/norm(cross(r,v));
S = cross(W,R);
Q = [R(1) S(1) W(1);R(2) S(2) W(2);R(3) S(3) W(3)];
p_ijk = Q*p_RSW;
p_dot_ijk = Q*p_dot_RSW;
r2 = r+p_ijk;
v2 = v+p_dot_ijk;

n = sqrt(mu/(a)^3);
omega = [0 0 n];
p_rel_noninert = p_RSW;
p_dot_rel = (p_dot_RSW'-cross(omega,p_RSW))*1000;

P = (2*pi)/n;
Y01 = [r v];
Y02 = [r2 v2];
tspan_sec = 0:60:3*P;
odeoptions = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T1,Y1] = ode45(@yprop, tspan_sec, Y01, odeoptions, mu);
[T2,Y2] = ode45(@yprop, tspan_sec, Y02, odeoptions, mu);
r1_ijk_km = Y1(:,1:3); % Target
v1_ijk_km = Y1(:,4:6);
r2_ijk_km = Y2(:,1:3); % Chaser
t = tspan_sec;
%% Modeling Target and Chaser Orbit Around the Earth
[X,Y,Z] = sphere;
rE = 6378;
figure(1)
hold on
xlabel('x Position (km)','fontsize', 12);
ylabel('y Position (km)','fontsize', 12);
zlabel('z Position (km)','fontsize', 12);
axis equal
view(-45,30); %to set viewing angle
%plot3(r2_ijk_km(:,1),r2_ijk_km(:,2),r2_ijk_km(:,3),'Color','none')
hSurface = surf(X*rE,Y*rE,Z*rE);
set(hSurface,'FaceColor',[0 0 1], 'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
xlim = [-9000 9000];
ylim = [-9000 9000];
zlim = [-9000 9000];
grid on 
filename = 'targetandchaser.gif';
p2 = plot3(r2_ijk_km(1,1),r2_ijk_km(1,2),r2_ijk_km(1,3),'g','LineWidth',2); %first iteration
m2 = scatter3(r2_ijk_km(1,1),r2_ijk_km(1,2),r2_ijk_km(1,3),'filled','g');
p1 = plot3(r1_ijk_km(1,1),r1_ijk_km(1,2),r1_ijk_km(1,3),'r','LineWidth',2); %first iteration
m1 = scatter3(r1_ijk_km(1,1),r1_ijk_km(1,2),r1_ijk_km(1,3),'filled','r');
legend('Earth','','Chaser','','Target')
for k = 1:length(t)
    p2.XData = r2_ijk_km(1:k,1); % line update
    p2.YData = r2_ijk_km(1:k,2); 
    p2.ZData = r2_ijk_km(1:k,3);
    m2.XData = r2_ijk_km(k,1); % point update
    m2.YData = r2_ijk_km(k,2);
    m2.ZData = r2_ijk_km(k,3);

    % for target debris trajectory
    p1.XData = r1_ijk_km(1:k,1); % line update
    p1.YData = r1_ijk_km(1:k,2);
    p1.ZData = r1_ijk_km(1:k,3);
    m1.XData = r1_ijk_km(k,1); % point update
    m1.YData = r1_ijk_km(k,2);
    m1.ZData = r1_ijk_km(k,3);
    
title(sprintf('Propogated Trajectory of Target and Chaser\nTime: %0.2f sec',t(k)),'fontsize', 12);
pause(0.1)

frame = getframe(gcf);
im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end
end

%% Relative Positon Derived from Inertial EOMs
r1_RSW_m = zeros(length(tspan_sec),3);
r2_RSW_m = r1_RSW_m;
p_rel = r1_RSW_m;

for i=1:length(tspan_sec)
 R = r1_ijk_km(i,:)/norm(r1_ijk_km(i,:));
 W = cross(r1_ijk_km(i,:),v1_ijk_km(i,:))/(norm(cross(r1_ijk_km(i,:),v1_ijk_km(i,:))));
 S = cross(W,R);
 Q = [R(1) S(1) W(1);R(2) S(2) W(2);R(3) S(3) W(3)];
 r1_RSW_m(i,1:3) = (Q'*r1_ijk_km(i,:)').*1000;
 r2_RSW_m(i,1:3) = (Q'*r2_ijk_km(i,:)').*1000;
 p_rel(i,1:3) = r2_RSW_m(i,:)-r1_RSW_m(i,:);
end

figure(3)
plot(p_rel(:,2),p_rel(:,1),'LineWidth',2) % plot in m
xlabel('y (m)','fontsize', 12);
ylabel('x (m)','fontsize', 12);
title('Relative Postions using Inertial EOMs','fontsize', 12);

%% Relative Position Derived from CW Targeting Equations
timespan = tspan_sec;
x0 = p_RSW(1);
y0 = p_RSW(2);
z0 = p_RSW(3);
x0dot = p_dot_rel(1)/1000;
y0dot = p_dot_rel(2)/1000;
z0dot = p_dot_rel(3)/1000;
for j = 1:length(timespan)
   t = timespan(j);
   x(j) = (4-3*cos(n*t))*x0 +(sin(n*t)/n)*x0dot + (2/n)*(1-cos(n*t))*y0dot;
   y(j) = 6*(sin(n*t)-n*t)*x0 + y0 + (2/n)*(cos(n*t)-1)*x0dot + (1/n)*(4*sin(n*t)-3*n*t)*y0dot;
   z(j) = cos(n*t)*z0 + (sin(n*t)/n)*z0dot;
pCW(j,:) = [x(j) y(j) z(j)].*1000;
end

figure(4)
hold on
plot(pCW(:,2),pCW(:,1),'LineWidth',2);
hold off
xlabel('y (m)','fontsize', 12);
ylabel('x (m)','fontsize', 12);
title('Relative Positons using CW Targeting','fontsize', 12 );

%% Error Between Targeting Methods
error = (p_rel-pCW)./1000;
figure(5)
subplot(3,1,1)
plot(timespan,error(:,1))
title('Error in Solutions via CW Equations')
subplot(3,1,2)
plot(timespan,error(:,2))
ylabel('True Relative Position - CW Derived Relative Position (m)')
subplot(3,1,3)
plot(timespan,error(:,3))
xlabel('Time (s)')
