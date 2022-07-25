%% This code simulates the deterministic bivirus model

clear variables
close all
clc

tic

%Declare global variables for use in the simulation
global D1 D2 B1 B2 n dfx A D

%% Obtain the infection matrix for virus 1
B1 = readmatrix('Total_RAW.csv'); B1 = B1(:,2:end); %Load adjacency matrix for virus 1

n = length(B1);   %Number of nodes in the network

D1 = 1.*eye(n);
D2 = 1.*eye(n);

% For convenience, we normalise B1 so that all rows sum to 2, and eliminate
% any edges with weight below a threshold epsilon.
P = 0.5*diag(B1*ones(n,1)); B1 = inv(P)*B1;

%This removes edges that are below epsilon value in weight
epsilon = 0.00005; C = double(B1>=epsilon);  
for i = 1:n
    for j = 1:n
        if C(i,j) == 0
            B1(i,j) = 0;
        end
    end
end

G = digraph(B1);  %Check the graph is still connected
[bins,binsizes] = conncomp(G);
if binsizes==n
    disp('G(B1) is strongly connected')
else
    disp('G(B1) is not strongly connected')
end

%Renormalise the matrix after deleting edges below threshold
P = 0.5.*diag(B1*ones(n,1)); B1 = inv(P)*B1;  


tfinal = 5000; tspan = [0 tfinal];  %Change 2nd term for simulation time

x0 = rand(n,1)'; A = B1; D = D1;
[t,x] = ode45(@sis,tspan,x0'); t = t'; x = x';   %Bivirus simulation
% [t,x] = ode45(@bivirus_fast,tspan,x0'); t = t'; x = x';   %Bivirus simulation for really slow time constant systems

x_bar = x(:,end);  %Set single virus endemic equilibrium


%% Construct B2
% The following code creates B2 which yields two stable boundary
% equilibria according to the PRE paper method

%Change k to change the vector e_i. We change the weights of edges incoming to node k.
k = 48; ei = zeros(n,1); ei(k) = 1;   
[a,b]= maxk(B1(k,:),2);   %Finds the two largest elements of row k in B1
z = zeros(n,1); z(b(1)) = 1; z(b(2)) = -x_bar(b(1))*z(b(1))/x_bar(b(2));  %Construct the vector z

eps = min([B1(k,b(1)),B1(k,b(2))])./(2*max(z));
Bp =  B1 + eps.*ei*z';
clear k a b

IX = eye(n) - diag(x_bar);
F = IX^(-2)-Bp;

[a,b] = eig(B1');   %Get left eigenvectors and eigenvalues of B1
[c,d] = eig(Bp');

[m,i] = max(real(diag(b)));  %Identify PF eigenvalue and index of B1
u = a(:,i);  u = u./(u'*x_bar);  %Get normalised left PF eigenvector of B1

[m,i] = max(real(diag(d)));  %Identify PF eigenvalue and index of Bp
v = c(:,i);  v = v./(v'*x_bar);  %Get normalised left PF eigenvector of Bp

ut = (u'*diag(x_bar)*inv(IX)*inv(F))';  %Tilde u
vt = (v'*diag(x_bar)*inv(IX)*inv(F))';  %Tilde v

[m,i] = max(ut./vt);   %Get largest difference of ut_i/vt_i

check5 = 0; check6 = 0;
counter = 0;
while (check5==0) || (check6 == 0)
    
    check5 = 0; check6 = 0;
    w = ones(n,1); w(i) = 0; w = w./sum(w);   %Uniform weights except at position i
    j = randsample(n,1,true,w);  %Select another index except i uniformly at random.
    
    k1 = vt(j)/vt(i);   %Smaller ratio  (note this is actually the inverse)
    k2 = ut(j)/ut(i);   %Larger ratio   (note this is actually the inverse)
    
    alpha = k2 + (k1-k2).*rand;    %Alpha should be between kj and ki;
    
    check1 = (alpha*ut(i))/ut(j) > 1;   %Check alpha works
    check2 = (alpha*vt(i))/vt(j) < 1;
    
%     j_ind = find(Bp(j,:));
%     i_ind = find(Bp(i,:));
%     
%     b_jk = randsample(j_ind,1);
%     b_ik = randsample(i_ind,1);

    [q,j_ind] = maxk(Bp(j,:),5);
    [q,i_ind] = maxk(Bp(i,:),5);
    
    b_jk = max(j_ind);
    b_ik = max(i_ind);
    
    beta = (0.5*Bp(j,b_jk)*rand)*x_bar(j);
    
    delta_b_jk = -beta/x_bar(j);
    delta_b_ik = (alpha*beta)/x_bar(i);
    
    delta_B = zeros(n);
    delta_B(j,b_jk) = delta_b_jk;
    delta_B(i,b_ik) = delta_b_ik;
    
    
    s = delta_B*x_bar;
    
    check3 = ut'*s < 1;
    check4 = vt'*s < 1;
    
    B2 = Bp+delta_B;
    
    delta_x = inv(F)*delta_B*x_bar;
    
    new_y_est = x_bar + delta_x;
    
    tspan = [0 1000];  %Change 2nd term for simulation time
    
    A = B2; D = D1;
    [t,x] = ode45(@sis,tspan,ones(n,1)); t = t'; x = x';   %Bivirus simulation
    
    new_y = x(:,end);
    q1 = max(real(eig((eye(n)-diag(new_y))*B1)));
    q2 = max(real(eig((eye(n)-diag(x_bar))*B2)));
%     check5 = max(eig((eye(n)-diag(new_y))*B1)) < 0.9999;
%     check6 = max(eig((eye(n)-diag(x_bar))*B2)) < 0.9999;
    check5 = q1 < 0.999999;
    check6 = q2 < 0.999999;
    
    min_q = min(min_q,q1+q2);
    
    clear x
    counter = counter + 1;  %Checks how many iterations of the algorithm has been run
end


IC = 3; %Initial conditions case number
switch IC
    case 1  %Used for PRE simulation
        
        x1 = 0.1+(0.3-0.1).*rand(n,1); x2 = 0.1+(0.6-0.1).*rand(n,1);
        %         x1 = 0.1+(0.6-0.1).*rand(n,1); x2 = 0.1+(0.3-0.1).*rand(n,1);
    case 2
        eta = 0.95;
        x1 = (0.5*eta).*ones(n,1); x2 = (1-eta).*ones(n,1);  %Bottom right
        
        %         eta = 0.95;
        %         x1 = (1-eta).*ones(n,1); x2 = (0.5*eta).*ones(n,1);  %Top left
        
    case 3
        
        b = rand(n,1);  a = rand(n,1); 
        x1 = a./(b+a); x2 = (0.1*b)./(a+b);  
        
        b = rand(n,1);  a = rand(n,1); 
        x1 = a./(b+a); x2 = (0.5*b)./(a+b); 
end
x0 = [x1', x2']';


%% Simulate bivirus dynamics

tic
tfinal = 150000;

tspan = [0 tfinal];  %Change 2nd term for simulation time


[t,x] = ode45(@bivirus,tspan,x0'); t = t'; x = x';   %Bivirus simulation
% [t,x] = ode45(@bivirus_fast,tspan,x0'); t = t'; x = x';   %Bivirus simulation for really slow time constant systems

x_final = x(:,end);
x1_final = x_final(1:n); x2_final = x_final(n+1:end);

%% Plot simulation output

fig = figure;
% p1 = plot(t,x(1:n,:),':b','LineWidth',1.5);
p1 = semilogx(t,x(1:n,:),':b','LineWidth',1.5);
hold on
% p2 = plot(t,x(n+1:end,:),'r','LineWidth',1.5);
p2 = semilogx(t,x(n+1:end,:),'r','LineWidth',1.5);
xlab = xlabel('Time,  t');
ylab = ylabel('Fraction of Infected');
leg = legend([p1(1) p2(1)],'Virus 1', 'Virus 2','box','off');
set(leg,'FontSize',14,'Location','Best')
set(xlab,'FontSize',12)
set(ylab,'FontSize',12)
axis([0 tfinal 0 1])

toc
