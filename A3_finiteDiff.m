%Solution to Ass 3 A2
N=12;

% Changed to the interval [-1, 1] from the notes 
h = 2/N;%2*pi/N;
h
x= (0:N)*h - 1;
x
ue = @(x) (1-x.^2).*exp(x);
f = @(x) (4.*x + 2).*exp(x);

%setup matrix A
A = zeros(N-1,N-1);

for i = 1:N-1
            if i+1 <= N-1
                A(i, i+1) = -16/(12*h^2);
            end
            
            A(i, i) = 30/(12*h^2) + 1;
            
            if i-1 <= N-1 && i-1 > 0
                A(i, i-1) = -16/(12*h^2) ;
            end
            
            if i +2 <= N-1
                   A(i, i+2) = 1/(12*h^2);
            end
            
            if i-2 <= N-1 && i-2 > 0             
                A(i, i-2) = 1/(12*h^2);
            end
end


A(1,1) = (10/(6*h^2) + 1); %2/(h^2) + 1;
A(1,2) = -1/ (2*h^2); %-1/h^2;
A(1,3) = -1/(3*h^2); %0;
A(1,4) = 1/(12*h^2);

A(N-1,N-1) = (10/(6*h^2) + 1); %2/(h^2) + 1;
A(N-1,N-2) = -1/ (2*h^2); %-1/h^2;
A(N-1,N-3) = -1/(3*h^2); %0;
A(N-1,N-4) = 1/(12*h^2);

% A(N-1,N-1) = 2/(h^2) + 1;
% A(N-1, N-2) = -1/h^2;
% A(N-1, N-3) = 0;

% A(N-2, N-2) = 1/(12*h^2);
A
X = x
X(1) = [];%x(1) - h
X(N) = [];%x(N+1) + h
% X(N-1) = []
X
% 
b = f(X);
% % b'
% 
U = A\b';

% X = [-1 X]
U = [0 U']
U = [U 0]
% 
plot(x, U, x, ue(x))
% 
% % ue(x)
maxerr = max(abs(U-ue(x)))




