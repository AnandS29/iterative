% Reference: page 279 of Basar's textbook, 'Dynamic noncooperative game theory'
x0 = [2;3];
Q1 = [0,0;
      0,1];
Q2 = [1,-1;
     -1, 1];
R1 = 1;
R2 = 1;
n = length(x0);
A = eye(2);
B1 = [1;
      0];
B2 = [0;
      1];
m1 = size(B1,2);
m2 = size(B2,2);
T = 20;
%%
K = cell(T,2);
Z = cell(T+1,2);
F = cell(T,1);
Z{T+1,1} = Q1;
Z{T+1,2} = Q2;
for_debug_use = cell(T,1);
for iter = 1:T
    t = T+1-iter;
    sol_K = [R1 + B1'*Z{t+1,1}*B1, B1'*Z{t+1,1}*B2; B2'*Z{t+1,2}*B1, R2 + B2'*Z{t+1,2}*B2]^(-1)*[B1'*Z{t+1,1}*A; B2'*Z{t+1,2}*A];
    for_debug_use{iter} = svd([R1 + B1'*Z{t+1,1}*B1, B1'*Z{t+1,1}*B2; B2'*Z{t+1,2}*B1, R2 + B2'*Z{t+1,2}*B2]);
    K{t,1} = sol_K(1,:);
    K{t,2} = sol_K(2,:);
    F{t} = A - B1 * K{t,1} - B2 * K{t,2};
    Z{t,1} = F{t}'*Z{t+1,1}*F{t} + K{t,1}'*R1*K{t,1} + Q1;
    Z{t,2} = F{t}'*Z{t+1,2}*F{t} + K{t,2}'*R2*K{t,2} + Q2;
end
%% roll out trajectory
x0 = [2;3];
x_traj = [x0, zeros(n, T)]; % the state trajectory
u1_traj = zeros(m1, T);
u2_traj = zeros(m2, T);
total_cost1 = 1/2*x0'*Q1*x0; % the total costs for the first agent
total_cost2 = 1/2*x0'*Q2*x0; % the total costs for the second agent
for t = 1:T
    u1_traj(:,t) = K{t,1}*x_traj(:,t);
    u2_traj(:,t) = K{t,2}*x_traj(:,t);
    x_traj(:,t+1) = A*x_traj(:,t) - B1*u1_traj(:,t) - B2*u2_traj(:,t);
    total_cost1 = total_cost1 + 1/2*x_traj(:,t+1)'*Q1*x_traj(:,t+1) + 1/2*u1_traj(:,t)'*R1*u1_traj(:,t);
    total_cost2 = total_cost2 + 1/2*x_traj(:,t+1)'*Q2*x_traj(:,t+1) + 1/2*u2_traj(:,t)'*R2*u2_traj(:,t);
end
total_cost1_close = total_cost1;
total_cost2_close = total_cost2;
fprintf('state trajectory:')
disp(x_traj)
