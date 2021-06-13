% set parameters
gamma = 20;

% set disturbances matrix
D = [0 0 0;
     1 0 0;
     0 0 0;
     0 1 0;
     0 0 0;
     0 0 1];

% Solve system of matrix equiations
cvx_begin sdp

    variable P(6,6) symmetric
    variable Y(1,6)
    P > 0;
    [A*P+P*A'+B*Y+Y'*B', D, P*C';
     D', -gamma*eye(3), zeros(3,6);
     C*P, zeros(6,3), -gamma*eye(6)] <= 0;
 
cvx_end

% Coeffs of controller
K = -Y*inv(P);
vpa(K)
