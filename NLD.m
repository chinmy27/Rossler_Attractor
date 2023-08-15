clear all
clc
close all
%% Parameters
A = [0.25 , 0.5 , 1 , 0.5];
B = [0.5 , 1 , 1.5 , 0.5];
a = 0.2;
b =0.2;
% for j = 1:length(A)
X_init = [2,2,0];
% a = A(j);
% b =B(j);
for c = 1:0.1:5


%% Fixed Points 
x1 = (c + (c^2 - 4*a*b)^(1/2))/2;
x2 = (c - (c^2 - 4*a*b)^(1/2))/2;
y1 = x1/-a;
y2 = x2/-a;
z2 = -y2;
z1 = -y1;
F1 = [x1;y1;z1];
F2 = [x2;y2;z2];
% 
% Jacobian
J1 = [0 -1 -1;1 a 0; z1 0 (x1-c)];
J2 = [0 -1 -1;1 a 0; z2 0 (x2-c)];
l1 = eig(J1);
l2 = eig(J2);

% Properties
T1 = trace(J1);
T2 = trace(J2);
D1 = det(J1);
D2 = det(J2);
T14 = T1^2 - 4*D1;
T24 = T2^2 - 4*D2;
%% RK

T = 0;
h = 0.01;

X(:,1) = X_init;
i = 1;

for T = 0:h:100
    k1 = h*Rossler(X(:,i),a,b,c);
    k2 = h*Rossler(X(:,i)+k1/2,a,b,c);
    k3 = h*Rossler(X(:,i)+k2/2,a,b,c);
    k4 = h*Rossler(X(:,i)+k3/2,a,b,c);
    X(:,i+1) = X(:,i) + 1/6*(k1+2*k2+2*k3+k4);

    i = i+1;
    
end
% subplot(2,2,j);
r = rem(c,2);
if c<2.5
plot3(X(1,:),X(2,:),X(3,:),'color','k');
else plot3(X(1,:),X(2,:),X(3,:),'color','r');
end
hold on
%plot3(F1(1),F1(2),F1(3),'ob',F2(1),F2(2),F2(3),'or');
% title('For Case ',j);
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
hold on
pause(0.01)
end

function Y = Rossler(X,a,b,c)
    x = X(1);
    y = X(2);
    z = X(3);
    Y(1) = -y - z;
    Y(2) = x + a*y;
    Y(3) = b + z*(x-c);
    Y = Y';
end

