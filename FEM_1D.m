% MATH 521, Homework 1, Problem A4
% FEM 1D implementation 
% The PDE is --- -u'' + u = f , u is 1-periodic in x

% exact solution 
u = @(x)  (sin(pi*x)).^2;

f = @(t) (-(2*pi^2*cos(pi*t).^2 - 2*pi^2*sin(pi*t)^2 ) +  (sin(pi*t)).^2);
% the grid 
N = 4;
h = 1/N;
x = (0:N)*h;
% x
% setup mass matrix M
M = zeros(N+1, N+1);
M00 = 2*h/3;
M11 = h/6;
M22 = h/6;
M = diag(M00*ones(1,N+1)) + diag(M11*ones(1,N),-1) + diag(M22*ones(1,N),1);
M(1, N+1) = M(1,1);
M(N+1,1) = M(N+1, N+1);
M

% setup stiffness matrix K
K = zeros(N+1, N+1);
K00 = 2/h;
K11 = -1/h;
K22 = -1/h;
K = diag(K00*ones(1,N+1)) + diag(K11*ones(1,N),-1) + diag(K22*ones(1,N),1);
K(1, N+1) = K(1,1);
K(N+1,1) = K(N+1, N+1); 

% setup the right hand side using quadrature
%f = zeros(N+1, 1);
%for i = 1:N+1
%    f(i) = u(x(i));
%end
F = zeros(N+1, 1);
for i = 2:N
     F(i) = (h/2)* (f(x(i) - (h/2)) + f(x(i) + (h/2)));
    h
    % for phi1 = (1 + xi)/2
    % g(-1/sqrt(3)) + g(1/sqrt(3))
    F1 = (h/2) *  ( (1-1/sqrt(3)) * f((h/2) * (-1/sqrt(3) + 1) + x(i-1)))  
    + (h/2) * ((1+1/sqrt(3)) * f((h/2) * (1/sqrt(3) + 1) + x(i-1)));
                  
    % for phi2 = (1 - xi)/2
    % g(-1/sqrt(3)) + g(1/sqrt(3))
    F2 = (h/2) *  ( (1+1/sqrt(3)) * f((h/2) * (-1/sqrt(3) + 1) + x(i)))  
    + (h/2) * ((1-1/sqrt(3)) * f((h/2) * (1/sqrt(3) + 1) + x(i)));
    
%     F(i) =  F1 + F2;
              
end
% % 
F(1) = (h/2) * (f(x(1) + (h/2)));
F(N+1) = (h/2)* (f(x(N+1) - (h/2)));

% 
% F(1) = (h/2) *  ( (1+1/sqrt(3)) * f((h/2) * (-1/sqrt(3) + 1) + x(1)))  
%     + (h/2) * ((1-1/sqrt(3)) * f((h/2) * (1/sqrt(3) + 1) + x(1)));
%   
% F(1) = 0;
% F(N+1) = 0 ; 
%     (h/2) *  ( (1-1/sqrt(3)) * f((h/2) * (-1/sqrt(3) + 1) + x(N)))  
%     + (h/2) * ((1+1/sqrt(3)) * f((h/2) * (1/sqrt(3) + 1) + x(N)));  

% F at the boundary -- you would have either of F1 or F2 at the boundary,
% but not both.

left = M + K;
right = F;

U = left\right;
% maxerr = max(abs(U-u(x)'));


x1 = linspace(0,1);
y1 = u(x1);%sin(4*pi*x1);
% 
figure(1)
%plot(x, h*fun1(x+h/2), x, fun1(x))
plot(x, U, x1, u(x1))
xlabel('x')
%ylabel('y')
legend('Approximation','Exact')

U = U(2:N)
x = x(2:N)
% convergence in max norm  O(h^2)
maxNormError = max(abs(U - u(x)'));
is = norm(U - u(x)', inf);
is2 = norm(U - u(x)', 'fro');

% calculate L2 error
L2_err = 0;

% for i=2:N
%     vec = [U(i), u(x(i))];
%     l2 = norm(vec, 'fro');
%     
%     diff = abs(U(i) - u(x(i)));
%     L2_err = L2_err + diff^2;
% end
% L2_err = sqrt(L2_err);
% L2_err


ue = u;
xp = x;
uh = U;
% x1 = 0:0.1:1;021 % for H1err we need the derivation of the exact solution
uprime = @(x) (2*pi*cos(pi*x)*sin*(pi*x));%str2func(['@(x)',char(diff(sym(ue)))]); 
% set the errors to zero, for summation
L2err = 0;
H1err = 0;
 
% h = diff(xp);
% duh = diff(uh)./h;
for j = 1:length(uh)
    %         L2err = L2err + sqrt(h(j))*(ue(xp(j))-uh(j))^2;
    L2err = L2err + (ue(xp(j))-uh(j))^2;
end

% for j = 1:length(duh)
%      %L2err = L2err + (ue(xp(j))-uh(j))^2;
%      H1err = H1err + ((uprime(xp(j)))-duh(j))^2;
% end
%  
% H1err = sqrt(mean(h))*sqrt(H1err + L2err);
L2err = sqrt(mean(h))*sqrt(L2err);
L2err

maxNormError

% for k = 1:length(xp)
%      uex(k) = ue(xp(k));
% end
%  
% LI0err = max(abs(uex-uh'));
% LI0err

% %y1 = 0:0.1:1
% % generate random grid 
% for i=1:11
%     % generate random number in the range 0,1
%     randomNum = rand()*2;
% 
%     % generate random number in the range [-h/3, h/3]
%     r = (2*h/3) * randomNum - (h/3);
%     r
%     
%     x1(i) = x1(i) + r;
% end
% x1
% y1 = zeros(11, 1);
% plot(x1, y1, 'o')
% % convergence in the HP1 norm
% 
% % convergence in the L2 norm
F