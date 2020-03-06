% Chebyshev spectral method
% The pde is -u'' + u = f
% boundary conditions = u(-1) = u(1) = 0

ue = @(x) (1-x.^2).*exp(x);
f = @(x) (4.*x + 2).*exp(x);

% number of intervals, order of the chebyshev interpolant
N = 8

% gets the chebyshev approximation of f
pf = chebfun(f, 'trunc', N);
format long
% get the cheyshev coefficients of the chebyshev approx - pf
bn = chebcoeffs(pf)

A = zeros(N+2, N);
b = zeros(N+2,1);
cn = 0 ;
p = 0;
for i = 1:N-2
        if i-1 == 0
           cn = 2;
        else
            cn = 1;
        end
        % the -u'' term
        p = 0;
        p = (i-1) + 2;
        while (p <= N-1)
            A(i, p+1) = -1/cn * (p *(p^2 - (i-1)^2) );
            p = p + 2;
        end
end

% the +u term
for i = 1:N
    A(i,i) = 1;
end

% boundary conditions 
for i = 1:N
    A(i,i)  = 1;
    A(N+1, i) = 1;
    if  mod(i-1,2) == 0
        A(N+2, i) = 1;
    else
        A(N+2, i) = -1;
    end
    
end

% construct b
for i = 1:N+2
    if i > N
        b(i,1) = 0;
    else
        b(i, 1) = bn(i);
    end
end

an = A\b;
an

% use chebfun to get the chebyshev coefficients of the known exact solution
uf = chebfun(ue, 'trunc', N);
format long
% get the cheyshev coefficients of the chebyshev approx - uf
un = chebcoeffs(uf)

% syms x
% cbt = chebyshevT(arr, x)

% chebyshev polynomials upto order 6 multiplied by the approximated
% chebyshev oefficients - an
fun = @(x) (an(1).*1 + an(2).*x + an(3).*(2.*x.^2 - 1) + an(4).*(4.*x.^3 - 3.*x) + an(5) * (8.*x.^4 - 8.*x.^2 + 1) + an(6) * (16.*x.^5 - 20.*x.^3 + 5.*x))

% plot the interpolant and exact function values, calculate error
x1 = -1:0.1:1;
plot(x1, fun(x1))
xlabel("x values on the grid")
ylabel("function values")
title("exact versus chebyshev interpolant of degree 6")
maxError = max(abs(ue(x1) - fun(x1)))
hold on
plot(x1, ue(x1))