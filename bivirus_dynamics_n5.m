%% This code simulates the deterministic bivirus model

clear variables
close all
clc

tic

%Declare global variables for use in the simulation
global D1 D2 B1 B2 n dfx A D

%% Parameter setup

n = 5;   %Number of nodes in the network

%Select the case parameters
recov = 1; %Diagonal recovery matrix case number
switch recov
    case 1
        D1 = 1.*eye(n);
        D2 = 1.*eye(n);
end

infect = 2; %B^i network infection matrix case number
switch infect
    case 1
        B1 = [1 0 0 0 1;1 1 0 0 0;0 1 1 0 0;0 0 1 1 0;0 0 0 1 1];
        
        % The following code creates B2 which yields two stable boundary
        % equilibria according to the PRE paper method
        Bp =  B1 + 0.5*[1 0 0 0 0]'*[-1 1 0 0 0];
        
        x_bar = 0.5.*ones(n,1);
        IX = eye(n) - diag(x_bar);
        F = IX^(-2)-Bp;
        
        [a,b] = eig(B1');   %Get left eigenvectors and eigenvalues of B1
        [c,d] = eig(Bp');
        
        [m,i] = max(diag(b));  %Identify PF eigenvalue and index of B1
        u = a(:,i);  u = u./(u'*x_bar);  %Get normalised left PF eigenvector of B1
        
        [m,i] = max(diag(d));  %Identify PF eigenvalue and index of Bp
        v = c(:,i);  v = v./(v'*x_bar);  %Get normalised left PF eigenvector of Bp
        
        ut = (u'*diag(x_bar)*inv(IX)*inv(F))';  %Tilde u
        vt = (v'*diag(x_bar)*inv(IX)*inv(F))';  %Tilde v
        
        [m,i] = max(ut./vt);   %Get largest difference of ut_i/vt_i
        
        check5 = 0; check6 = 0;
        
        while (check5==0) || (check6 == 0)
            
            check5 = 0; check6 = 0;
            w = ones(n,1); w(i) = 0; w = w./sum(w);   %Uniform weights except at position i
            j = randsample(n,1,true,w);  %Select another index except i uniformly at random.
            
            k1 = vt(j)/vt(i);   %Smaller ratio  (note this is actually the inverse)
            k2 = ut(j)/ut(i);   %Larger ratio   (note this is actually the inverse)
            
            alpha = k2 + (k1-k2).*rand;    %Alpha should be between kj and ki;
            
            check1 = (alpha*ut(i))/ut(j) > 1;   %Check alpha works
            check2 = (alpha*vt(i))/vt(j) < 1;
            
            j_ind = find(Bp(j,:));
            i_ind = find(Bp(i,:));
            
            b_jk = randsample(j_ind,1);
            b_ik = randsample(i_ind,1);
            
            
            beta = (0.1*Bp(j,b_jk)*rand)*x_bar(j);
            
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
            
            tspan = [0 10000];  %Change 2nd term for simulation time
            
            A = B2; D = D1;
            [t,x] = ode45(@sis,tspan,ones(n,1)); t = t'; x = x';   %Bivirus simulation
            
            new_y = x(:,end);
            q1 = max(eig((eye(n)-diag(new_y))*B1));
            q2 = max(eig((eye(n)-diag(x_bar))*B2));
            check5 = max(eig((eye(n)-diag(new_y))*B1)) < 0.9999;
            check6 = max(eig((eye(n)-diag(x_bar))*B2)) < 0.9999;
            
            clear x
        end
        
        
        
    case 2   %For PRE paper, n = 5 simulation
        B1 = [1 0 0 0 1;1 1 0 0 0;0 1 1 0 0;0 0 1 1 0;0 0 0 1 1];
        B2 = [0.500000000000000,0.500000000000000,0,0,1;0.934952590675814,1,0,0,0;0,1.08955953684245,1,0,0;0,0,1,1,0;0,0,0,1,1];
        
end

IC = 5; %Initial conditions case number
switch IC
    case 1
        b = rand(n,1);  a = rand(n,1); c = rand(n,1);
        x1 = a./(b+a+c); x2 = b./(a+b+c);
    case 2
        x1 = rand(n,1);  x2 = 0.00.*ones(n,1);
        %         x2 = rand(n,1);  x1 = 0.00.*ones(n,1);

    case 3
        x2 = rand(n,1); x2(1:3) = 0;  x1 = 0.*ones(n,1);
    case 4
        eta = 0.1;
        x1 = (0.55*eta).*ones(n,1); x2 = (0.5*eta).*ones(n,1);
        
%                 x2 = (0.55*eta).*ones(n,1); x1 = (0.3*eta).*ones(n,1);
    case 5  %Used for PRE simulation

        x1 = 0.1+(0.3-0.1).*rand(n,1); x2 = 0.1+(0.6-0.1).*rand(n,1);
%         x1 = 0.1+(0.6-0.1).*rand(n,1); x2 = 0.1+(0.3-0.1).*rand(n,1);
    case 6
        eta = 0.95;
        x1 = (0.5*eta).*ones(n,1); x2 = (1-eta).*ones(n,1);  %Bottom right
        
%         eta = 0.95;
%         x1 = (1-eta).*ones(n,1); x2 = (0.5*eta).*ones(n,1);  %Top left
end
x0 = [x1', x2']';


%% Simulate bivirus dynamics

s1 = max(real(eig(-D1+B1)));   disp(['s1 = ',num2str(s1)])
s2 = max(real(eig(-D2+B2)));   disp(['s2 = ',num2str(s2)])
p1 = max(eig(inv(D1)*B1)); disp(['rho[(D^1)^{-1}B^1] = ',num2str(p1)])
p2 = max(eig(inv(D2)*B2)); disp(['rho[(D^2)^{-1}B^2] = ',num2str(p2)])


tlength = 1500;
tfinal = 50000;

tspan = [0 tfinal];  %Change 2nd term for simulation time

% tspan = 0:0.2:tfinal;
% x = zeros(2*n,length(tspan));

[t,x] = ode45(@bivirus,tspan,x0'); t = t'; x = x';   %Bivirus simulation
% [t,x] = ode45(@bivirus_fast,tspan,x0'); t = t'; x = x';   %Bivirus simulation for really slow time constant systems

x_final = x(:,end);
x1_final = x_final(1:n); x2_final = x_final(n+1:end);

dfx = jacob(x(1:n,end),x(n+1:end,end));
dfx_eig = eig(dfx);  %Eigenvalues of the Jacobian at steady state
disp(['Largest real part of Jacobian eigenvalue: ',num2str(max(real(eig(dfx))))])

p3 = max(eig(inv(D1)*(eye(n)-diag(x2_final))*B1)); disp(['rho[(D^1)^{-1}(I-Z^2)B^1] = ',num2str(p3)])
p4 = max(eig(inv(D2)*(eye(n)-diag(x1_final))*B2)); disp(['rho[(D^2)^{-2}(I-Z^1)B^2] = ',num2str(p4)])

%% Plot simulation output

figure
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
