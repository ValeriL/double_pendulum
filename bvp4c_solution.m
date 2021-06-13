% clear all;

syms t

% unknown parameters
p1 = 1;
p2 = 1;
p3 = 1;
p4 = 1;
p5 = 1;
p6 = 1;

p = [p1 p2 p3 p4 p5 p6];
w = inv_T;

% problem solution
xmesh = linspace(0,swing_up_time,30);
solinit = bvpinit(xmesh,@guess,p);
solution = bvp4c(@system, @boundaries, solinit);

% input function parameters
P = solution.parameters;

xint = linspace(0,swing_up_time);
Sxint = deval(solution,xint);

plot(xint,Sxint)
xlabel('Time (seconds)')
ylabel('y')
legend('x','v', '{\theta_1}', '{w_1}', '{\theta_2}', '{w_2}')

function u = input(t,p)

    w = inv_T;
    p1 = p(1)*cos(2*w*t);
    p2 = p(2)*cos(3*w*t);
    p3 = p(3)*cos(4*w*t);
    p4 = p(4)*cos(5*w*t);
    p5 = p(5)*cos(6*w*t);
    p6 = p(6)*cos(7*w*t);
    a0 = -p(1)-p(3)-p(5);
    a1 = -p(2)-p(4)-p(6);
    u = a0+a1*cos(w*t)+p1+p2+p3+p4+p5+p6;
end

function dxdt = system(t,x,p)

    % system descrption
    l1 = 1;
    l2 = 1;
    m1 = 0.9;
    m2 = 1.1;
    m = 3; 
    g = 9.8;
    d1 = 0.2;
    d2 = 0.01;
    d3 = 0.01;
    
    % input
    u = input(t,p);
    
    % system ode
    dxdt(1)=x(2);
    dxdt(2)=(m1*(l2*m2*x(6)^2*sin(2*x(3)-x(5))-g*m2*sin(2*x(3))-g*m1*sin(2*x(3))+2*l1*m1*x(4)^2*sin(x(3))+2*l1*m2*x(4)^2*sin(x(3))+l2*m2*x(6)^2*sin(x(5))))/(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5)))+(u*(2*m1+m2-m2*cos(2*x(3)-2*x(5))))/(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5)))-(d2*l2*m2*x(4)*cos(x(3)-2*x(5))+2*d1*l1*l2*m1*x(2)+d1*l1*l2*m2*x(2)+d3*l1*m1*x(6)*cos(2*x(3)-x(5))+d3*l1*m2*x(6)*cos(2*x(3)-x(5))-2*d2*l2*m1*x(4)*cos(x(3))-d2*l2*m2*x(4)*cos(x(3))-d3*l1*m1*x(6)*cos(x(5))-d3*l1*m2*x(6)*cos(x(5))-d1*l1*l2*m2*x(2)*cos(2*x(3)-2*x(5)))/(l1*l2*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))));
    dxdt(3)=x(4);
    dxdt(4)=-(l1*m1^2*x(4)^2*sin(2*x(3))-2*g*m*m1*sin(x(3))-g*m*m2*sin(x(3))-2*g*m1*m2*sin(x(3))-2*g*m1^2*sin(x(3))-g*m*m2*sin(x(3)-2*x(5))+l2*m1*m2*x(6)^2*sin(x(3)+x(5))+2*l2*m*m2*x(6)^2*sin(x(3)-x(5))+l2*m1*m2*x(6)^2*sin(x(3)-x(5))+l1*m1*m2*x(4)^2*sin(2*x(3))+l1*m*m2*x(4)^2*sin(2*x(3)-2*x(5)))/(l1*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))))-(2*d2*l2*m*x(4)+2*d2*l2*m1*x(4)+d2*l2*m2*x(4)-2*d3*l1*m*x(6)*cos(x(3)-x(5))-d3*l1*m1*x(6)*cos(x(3)-x(5))-d3*l1*m2*x(6)*cos(x(3)-x(5))-d2*l2*m2*x(4)*cos(2*x(5))+d3*l1*m1*x(6)*cos(x(3)+x(5))+d3*l1*m2*x(6)*cos(x(3)+x(5))-2*d1*l1*l2*m1*x(2)*cos(x(3))-d1*l1*l2*m2*x(2)*cos(x(3))+d1*l1*l2*m2*x(2)*cos(x(3)-2*x(5)))/(l1^2*l2*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))))-(u*(2*m1*cos(x(3))+m2*cos(x(3))-m2*cos(x(3)-2*x(5))))/(l1*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))));
    dxdt(5)=x(6);
    dxdt(6)=(m*(g*m1*sin(x(5))+g*m2*sin(x(5))-g*m1*sin(2*x(3)-x(5))-g*m2*sin(2*x(3)-x(5))+2*l1*m1*x(4)^2*sin(x(3)-x(5))+2*l1*m2*x(4)^2*sin(x(3)-x(5))+l2*m2*x(6)^2*sin(2*x(3)-2*x(5))))/(l2*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))))-(d3*l1*m1^2*x(6)+d3*l1*m2^2*x(6)+d2*l2*m2^2*x(4)*cos(x(3)+x(5))-d2*l2*m2^2*x(4)*cos(x(3)-x(5))+2*d3*l1*m*m1*x(6)+2*d3*l1*m*m2*x(6)+2*d3*l1*m1*m2*x(6)-d3*l1*m1^2*x(6)*cos(2*x(3))-d3*l1*m2^2*x(6)*cos(2*x(3))+d2*l2*m1*m2*x(4)*cos(x(3)+x(5))+d1*l1*l2*m2^2*x(2)*cos(2*x(3)-x(5))-2*d2*l2*m*m2*x(4)*cos(x(3)-x(5))-d2*l2*m1*m2*x(4)*cos(x(3)-x(5))-d1*l1*l2*m2^2*x(2)*cos(x(5))-2*d3*l1*m1*m2*x(6)*cos(2*x(3))+d1*l1*l2*m1*m2*x(2)*cos(2*x(3)-x(5))-d1*l1*l2*m1*m2*x(2)*cos(x(5)))/(l1*l2^2*m2*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))))+(u*(cos(2*x(3)-x(5))-cos(x(5)))*(m1+m2))/(l2*(2*m*m1+m*m2+m1*m2-m1^2*cos(2*x(3))+m1^2-m1*m2*cos(2*x(3))-m*m2*cos(2*x(3)-2*x(5))));
    dxdt=dxdt';
end

% Начальная и конечная точки
function res = boundaries(ya,yb,p)
    res=[ya(1) yb(1) ya(2) yb(2) ya(3)-pi yb(3) ya(4) yb(4) ya(5)-pi yb(5) ya(6) yb(6)];
end

% Предположительные траектории
function g = guess(t)   

    w = inv_T;
    g = [0 0 pi-w*t 0 pi-w*t 0];
end


function w = inv_T
    w = pi/swing_up_time;
end

% desired swing-up time
function T = swing_up_time
    T = 3.3;
end
