
syms theta1 theta2 w1 w2 input f1 f2 f3 qpr l1 l2 m m1 m2 g d1 d2 d3

%% Non-linerized model

% M*diff(diff(y)) = f - d*diff(y) + u + w

M = [m+m1+m2 l1*(m1+m2)*cos(theta1) m2*l2*cos(theta2);
    l1*(m1+m2)*cos(theta1) l1^2*(m1+m2) l1*l2*m2*cos(theta1-theta2);
    l2*m2*cos(theta2) l1*l2*m2*cos(theta1-theta2) l2^2*m2];

F = [l1*(m1+m2)*w1^2*sin(theta1)+m2*l2*w2^2*sin(theta2);
    -l1*l2*m2*w2^2*sin(theta1-theta2)+g*(m1+m2)*l1*sin(theta1);
    l1*l2*m2*w1^2*sin(theta1-theta2)+g*l2*m2*sin(theta2)];

% input
input = [input; 0; 0];

% disturbances
f = [f1; f2; f3];

% damping
d = [d1*qpr; d2*w1; d3*w2];

Yprpr = simplify(M^(-1)*F);
dprpr = simplify(M^(-1)*d);
uprpr = simplify(M^(-1)*input);
fprpr = simplify(M^(-1)*f);

% Second-order system of equations
Final = Yprpr - dprpr + uprpr + fprpr;

%% Linearized model

% Variables values
l1 = 1;
l2 = 1;
m1 = 0.9;
m2 = 1.1;
m = 3; 
g = 9.8;
d1 = 0.2;
d2 = 0.01;
d3 = 0.01;

% State-space model
A = [0 1 0 0 0 0;
    0 -d1/(m*m1) -g*(m2+m1)/m d2/(l1*m) 0 0;
    0 0 0 1 0 0;
    0 d1/(l1*m) g*(m*m1+m*m2+m1*m2+m1^2)/(m*m1*l1) -d2*(m+m1)/(l1^2*m*m1) -g*m2/(m1*l1) d3/(l1*l2*m1);
    0 0 0 0 0 1;
    0 0 -g*(m1+m2)/(l2*m1) d2/(l1*l2*m1) g*(m1+m2)/(l2*m1) -d3*(m1+m2)/(l2^2*m2*m1)];

B = [0;
    1/m;
    0;
    -1/(l1*m);
    0;
    0];

C = eye(6);

% Check controllability
disp("Is controllable?")
if rank(ctrb(A,B)) == 6
    disp("Yes")
else
    disp("No")
end

%% LQR

% Larger coeff - larger penalized
Q = [150 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 300 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 300 0;
    0 0 0 0 0 1];

% A coefficient of u, penalizes input
R = 0.1;

% Coeffs of controller
K = lqr(A,B,Q,R);

%% Модальное управление

% desired time of transient period
tp_star = 11.83;
tp = 5;

w0 = tp_star/tp;

% Desired roots
poly_coeffs = [1 3.86*w0 7.46*w0^2 9.13*w0^3 7.46*w0^4 3.86*w0^5 w0^6];

% Pole placement
poles = roots(poly_coeffs);
K = place(A,B,poles);

% Check eigvalues 
A_new = A-B*K;
e = eig(A);
e_new = eig(A_new)

% Check response
% sys_new = ss(A_new,B,C,D);
% step(sys_new)

% Plot pole-zero map
D = [0];
pzmap(ss(A,B,C,D),ss((A-B*K),B,C,D),'r')
legend('initial','with regulator')
