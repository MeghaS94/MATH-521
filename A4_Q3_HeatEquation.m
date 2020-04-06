% MATH 521, Homework 4, Problem A3
% FEM Time dependent - heat equation
% The PDE is --- du/dt =  u'' , u(0) = u(1) = 0
% I also need an initial condition ..

% exact solution 
ue = @(x,t) cos(x)*exp(-t)

ns = []
max_errs = []
l2_errs = []

for t = 0:8
    % the grid
    N = 16 %*(2^t);
    % Changed to the interval [0, 1]
    h = 1/N;
    x = (0:N)*h*pi - pi/2;
    h =h*pi;
%     x
%     ue(x,1)

    i = t;
    r = (rand(N+1, 1)*2*(h/3^i) - (h/3^i))';
%     x = x + r;
%     x_r = [-pi/2 x]
    x_r(1) = -pi/2;
    x_r(N+2) = pi/2; 
    % plot(x, zeros(9, 1), 'o')

    % setup the mass matrix 
    M = zeros(N+1, N+1);
    % setup stiffness matrix K
    K = zeros(N+1, N+1);

    indx = 2;
    for i =1:N+1
        if indx <N+1
             x1 = x_r(indx-1);
             if indx + 1 <= N+1
                x2 = x_r(indx);
                 x3 = x_r(indx +1);
             end
             if indx <=N+1
                 x2 = x_r(indx);
             end
        end

        for j = 1:N+1
            if i == j 
                if i == 1 
                    M(i,j) = 2*(x_r(2) - x_r(1)) / 3; 
                    M(i, i+1) = (x_r(2) - x_r(1))/6;
                    K(i,i+1) = -1/(x_r(2) - x_r(1));
                     K(i,j) = 2 / (x_r(2) - x_r(1));
                elseif i == N+1
                    M(i,j) = 2*(x_r(N+1) - x_r(N)) / 3; 
                    M(i, i-1) = (x_r(N+1) - x_r(N))/6;
                    K(i,i-1) = -1/(x_r(N+1) - x_r(N));   
                     K(i,j) = 2 / (x_r(N+1) - x_r(N));
                else
                    M(i,j) = (x2-x1) / 3 + (x3-x2)/3; 
                    M (i, j-1) = (x2-x1)/6;
                    M (i,j+1) = (x3-x2)/6;

                    K(i, j-1) = -1/(x2-x1);
                    K(i,j+1) = -1/(x3-x2);
                    K(i,j) = 1 / (x2-x1) + 1/(x3-x2);
                end
            end
        end
        indx = indx +1;
    end



    % M2 = M;
    % mass lumping
    for i = 1:N+1
        for j = 1:N+1
            if i == j 
                if i ==1
                    M(1,1) = M(1,1) + 2*M(1,2);
                    M(1,2) = 0;
                elseif i == N+1
                    M(N+1, N+1) = M(N+1, N+1) + 2*M(N+1, N);
                    M(N+1, N) = 0;
                else
                    M(i,i) = M(i,i) + M(i, i-1) + M(i, i+1);
                    M(i, i-1) = 0;
                    M(i, i+1) = 0;
                end
            end
        end
    end
%     M
    % M(1,1) = M(1,1)/2;
    % M(N+1,N+1) = M(N+1,N+1)/2;

    K(1,N+1) = K(1,2);
    K(N+1,1) = K(N+1,N);
%     K

    % finite difference matrix 
    e = ones(N+1,1);
    A = spdiags([e -(2)*e e], -1:1, N+1, N+1);
    % boundary conditions 
    % A(1,2) = 2;
    % A(N+1, N) = 2;

    A(1,N+1) = 1;
    A(N+1, 1) = 1;

    % sanity checks
    a = full(A)/h^2;

    % FEM matrix
    A2 = inv(M)*K;
    % A2(1,1) = A2(1,1)/2;
    % A2(N+1,N+1) = A2(N+1,N+1)/2;

    % initial conditions
    U0 = cos(x);
    U0 = U0';
    % or
    % U0 = rand(N+1, 1); %ones(N+1,1)*rand(); 

    k = 0.0005;
    m = floor(1/k)+1;
    k = 1/m;
    t_current = 0;
    k1 = k/(h^2);
    C = 1/18;
    k2 = 0.04*h^2;
%     k2
    % plot(x,U0)
    % forward euler on the heat equation 
    for j = 2:m % for time steps
        U = U0;   
        U = U - k*A2*U;
        t_current = (j-1)*k;

    %     The FD approximation
    %     U = U + k1*A*U;
    %     t_current = (j-1)*k;


    %         t_current = (j-1)*k2;
        for i = 1:N+1 % for space steps
            if i == 1 || i == N+1   
                U(i) = 0; 
            elseif i == N+1
    %             U(i) = U(N);
            else
    %              U(i) = U(i)+k*(U(i+1)-2*U(i)+U(i-1));
            end
        end

        U0 = U;
%         plot(x,U0, 'b')
%         hold on;- '*')
%         ylim([-2, 2])
% 
%          pause(0.05)   
%          clf
    end
%     t_current
    % Error at time T - wrt grid size
    % max norm error

    % plot(x, ue(x, t_current), '*')
    % ylim([-2, 2])
    % hold on
    % plot(x, U, 'b')

    max_norm_error = max(abs(U - ue(x, t_current)' ))
    % ue(x, t_current)
    % U'
     
     uh = U;

     % set the errors to zero, for summation
     L2err = 0;

     for j = 2:length(uh)-1
        L2err = L2err + (ue(j, t_current)-uh(j)^2);
    end

    L2err = sqrt(mean(h))*sqrt(L2err)
%     ns = [ns, N];
    ns = [ns t];
    max_errs = [max_errs, max_norm_error];
    l2_errs = [l2_errs, L2err];
end

ns
max_errs
l2_errs
% 
plot(log(ns), max_errs) %legend('L2', 'Max Norm')
title("With decreasing values of h on a non-uniform grid")
xlabel('log(N)')
ylabel('max norm error')


% plot(ones(t+1,1)'*t - ns, max_errs) %legend('L2', 'Max Norm')
% title("With decreasing perturbation on a non-uniform grid")
% xlabel('degree of perturbation')
% ylabel('max norm error')

