% Reference: page 269 of Basar's textbook, 'Dynamic noncooperative game theory'
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
M = cell(T,2);
M{T,1} = Q1;
M{T,2} = Q2;
lambda = cell(T,1);
for iter = 1:T-1
    t = T-iter;
    lambda{t} = eye(n) + B1*R1^(-1)*B1'*M{t+1,1} + B2*R2^(-1)*B2'*M{t+1,2};
    M{t,1} = Q1 + A'*M{t+1,1}*lambda{t}^(-1)*A;
    M{t,2} = Q2 + A'*M{t+1,2}*lambda{t}^(-1)*A;
end
x_traj = [x0,zeros(n,T-1)]; % state trajectory
u1_traj = zeros(m1,T);
u2_traj = zeros(m2,T);
total_cost1 = 1/2*x0'*Q1*x0; % the total costs for the first agent
total_cost2 = 1/2*x0'*Q2*x0; % the total costs for the second agent
for t = 1:T-1
    u1_traj(:,t) = -R1^(-1)*B1'*M{t+1,1}*lambda{t}^(-1)*A*x_traj(:,t);
    u2_traj(:,t) = -R2^(-1)*B2'*M{t+1,2}*lambda{t}^(-1)*A*x_traj(:,t);
    x_traj(:,t+1) = lambda{t}^(-1)*(A*x_traj(:,t));
    total_cost1 = total_cost1 + 1/2*x_traj(:,t+1)'*Q1*x_traj(:,t+1) + 1/2*u1_traj(:,t)'*R1*u1_traj(:,t);
    total_cost2 = total_cost2 + 1/2*x_traj(:,t+1)'*Q2*x_traj(:,t+1) + 1/2*u2_traj(:,t)'*R2*u2_traj(:,t);
end
total_cost1_open = total_cost1;
total_cost2_open = total_cost2;

fprintf('state trajectory:')
disp(x_traj)