% A5 - FEM solver for a non-linear problem
% using (vector)newtons iterations

% The problem :
% -u'' + u^3 = f(x)
% u is periodic

ue = @(x)sin(pi*x)
f = @(x) pi^2*(sin(pi*x))  + sin(pi*x).^3

% checking if f has the correct value
% x = -5:0.1:5 
% plot(x,f(x))

% the grid 
N = 2000 ;
h = 1/N;
x = (0:N)*h;
% x

% staring guess for the solution U
U = ones(N+1, 1);

% stiffness martrix
K = zeros(N+1, N+1);
K00 = 2/h;
K11 = -1/h;
K22 = -1/h;
K = diag(K00*ones(1,N+1)) + diag(K11*ones(1,N),-1) + diag(K22*ones(1,N),1);
% TO DO - do I need boundary conditions here?
K(1, N+1) = K(1,1);
K(N+1,1) = K(N+1, N+1); 
% K

% evaluate non-linear term
I = zeros(N+1, 1);
for i= 2:N
    U0 = U(i-1);
    U1 = U(i);
    U2 = U(i+1);
    I(i) = h/2 * ( ((U0 + U1)/2)^3  + ((U1 + U2)/2)^3 );
end
% I
% TO DO - what do I do at the boundaries.....

% evaluate right hand side with quadrature
F = zeros(N+1, 1);
for i = 2:N
%     temp = f(x(i) - h/2)
     F(i) = (h/2)* (f(x(i) - (h/2)) + f(x(i) + (h/2)));
end
% 
F(1) = (h/2) * (f(x(1) + (h/2)));
F(N+1) = (h/2)* (f(x(N+1) - (h/2)));
% TO DO - (check if boundary conditions are consistent) F at the boundary -- you would have either of F1 or F2 at the boundary,
% but not both.
% F

% Non-linear function N(U)
NL = @(U) K*U + I - F;

% residual tolerance
r_tol = 1e-10;

% initial residual 
residual = NL(U);

% flag to detect when newton's iterations are diverging
fail = 0;
% iteration number
it = 0;

% newton iterates
while  abs(max(residual)) > r_tol && fail == 0
%     it
    Un = -jac(K,N,h,U)\NL(U);
    U = U + Un;
%     U
    residual = NL(U);
    it = it + 1;
end

% fem = K\F
plot(x, U, x, ue(x))
% U(1) =0;
% U(N) =0;
% F
% ue(x)
% plot(ue(x))
uevec = ue(x)';
uerr = max(abs(U-uevec));
fprintf('N is - %d and the error is - %d \n',N, uerr);


function J = jac(K,N,h,U) 
    % Jacobian matrix 
    J = zeros(N+1, N+1);
    J = J + K;

    % derivative of non-linear term 
    D = zeros(N+1, N+1);
    for i=2:N
         U0 = U(i-1);
        U1 = U(i);
        U2 = U(i+1);
        D(i,i) = (3*h/4)* ( ((U0+U1)/2)^2 + ((U1+U2)/2)^2  )  ;
        D(i, i-1) = (3*h/4)* ( ((U0+U1)/2)^2 );
        D(i, i+1) = (3*h/4)* ( ((U1+U2)/2)^2) ;
    end
    % TO DO - What do I do for the first and last row - boundaries
    D(1,1) = (3*h/4)* ( ((U(1)+U(2))/2)^2  ) ;
    D(1,2) = (3*h/4)* ( ((U(1)+U(2))/2)^2);
    D(N+1, N+1) = (3*h/4)* ( ((U(N+1)+U(N-1))/2)^2  ) ;
    D(N+1, N) = (3*h/4)* ( ((U(N+1)+U(N-1))/2)^2  );
    % D
    J = J + D;
end

