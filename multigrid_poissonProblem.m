% Assignment 5 Question A4
% Multi-grid method for the finite difference discretization of
% the 2D poisson operator

% Test problem in 2D using a uniform grid
N= 128;
pow = 7; 
Nt = N*N;
Nt_temp = (N+1)*(N+1);
h=1/N;

x=(0:N)*h;
y=(0:N)*h;
[X, Y] = meshgrid(x,y);
    
ct = cos(2*pi*X);
st = sin(2*pi*X);
uexact = exp(ct+Y).*Y.*(1-Y);
f = exp(ct+Y).*(-Y.^2-3*Y + ...
    4*pi*pi*(st.*st-ct).*(Y-Y.^2));

f = f';
fvec = f(:);
uexact = uexact';
uevec = uexact(:);

% surf(X,Y, uexact)
% hold on;
% Finite difference discretization of the poission problem
A = finiteDiffGrid(N+1);
size(A)
% For debugging
fullA = full(A);
D = spdiags(diag(full(A)), 0, Nt_temp, Nt_temp);
% full(D)

u = A\fvec;
ufd = u;
% This is the approximate solution on the 2D grid
ugrid = reshape(u, N+1, N+1);
uguess = zeros(Nt_temp, 1);



% 2-grid correction scheme
umgvec = MG(-A, -D, uguess, fvec, h, pow, ufd, 10000, X, Y);
% size(umg)
% max(abs(ujacobi - umg))
umg = reshape(umgvec, N+1, N+1);

% colormap([1 0 0;0 0 1]) %red and blue
% surf(X, Y, ugrid, 'FaceColor','g' )
% hold on;
% surf(X, Y, -umg, 'FaceColor','b' )


fprintf("MG error %d", max(abs(u+umgvec)));
% surf(X,Y, ugrid+umg)
% size(A)

% functions defined below
% ----------------------------------------------------------------- 
% function for jacobi iterations 
function u = jacobi(A, D, u, f, numIter)
    omega = 0;
    for i = 1:numIter
        u = omega*u + (1-omega) * inv(D)* (f - (A-D)*u);
    end
end

% Here u comes in as a 1D vector
function u = MG(A, D, u0, f, h, n, ufd, Njac, X, Y)
    % step 1 : Pre-smoothing/relaxing (apply jacobi iterations)
     ujacobi = jacobi(A, D, u0, f, Njac);
     ugrid = reshape(ujacobi, 2^n+1, 2^n+1)';
     
    % step 2 : Calculate the fine grid residual 
    r = A*ujacobi - f; %(TO DO : Check if the residual sign needs to be flipped)
    r_fine_grid = reshape(r, 2^n+1, 2^n+1)';
    fprintf("Error and residual at grid level 1-------   %d, %d \n", max(abs(ufd+ujacobi )), max(abs(r)))  
    
    E1 = reshape((ufd+ujacobi), 2^n+1, 2^n+1);
%     E1 = reshape(r, 2^n+1, 2^n+1);
    surf(X, Y, E1', 'FaceColor','r')
    hold on;
    
    % step 3 : restrict fine grid residual to coarse grid
    r_coarse =  restriction(r_fine_grid, n-1);
    % temp array becase I need the row major ordering in the vector
    r_coarse_temp = r_coarse';
    r_coarse_vec = r_coarse_temp(:);
   
        % step 4 : solve the residual equation on the coarse grid 
        % I need A on the restricted grid
        A_coarse = finiteDiffGrid(2^(n-1)+1);
%         ac = size(A_coarse)
        % get coarse grid error
        e_coarse = A_coarse\r_coarse_vec; % --------> here I do direct solve, but you could do jacobi too.
        e_coarse_grid = reshape(e_coarse, 2^(n-1)+1, 2^(n-1)+1)';
        
        e_c_jacobi = jacobi(A_coarse, spdiags(diag(full(A_coarse)), 0, size(A_coarse,1),  size(A_coarse,1)), zeros(size(r_coarse_vec,1),1) , r_coarse_vec, Njac);
        e_c_grid = reshape(e_c_jacobi, 2^(n-1)+1, 2^(n-1)+1)';
          H = 1/(2^(n-1)+1);
                x=(0:2^(n-1))*H;
                y=(0:2^(n-1))*H;
                size(x)
                [x, y] = meshgrid(x,y);
        %         size(e_c_grid)
                surf(x, y,  e_c_grid , 'FaceColor','g')

        fprintf("Error and residual at grid level 2 ------- %d, %d \n", max(abs(e_c_jacobi )), max(abs(A_coarse*e_c_jacobi-r_coarse_vec)))  

                % calculate residual - prepare to move to grid3
                r1 = A_coarse*e_c_jacobi - r_coarse_vec;
%                 %         E1 = reshape(r, 2^n+1, 2^n+1);
%                 H = 1/(2^(n-1)+1);
%                 x=(0:2^(n-1))*H;
%                 y=(0:2^(n-1))*H;
%                 size(x)
%                 [x, y] = meshgrid(x,y);
%         %         size(e_c_grid)
%                 surf(x, y, reshape(r1, 2^(n-1)+1, 2^(n-1)+1) )
                
                r1_fine_grid = reshape(r1, 2^(n-1)+1, 2^(n-1)+1);

                % grid 3
                r1_coarse = restriction(r1_fine_grid, n-2);
                r1_coarse_temp = r1_coarse';
                r1_coarse_vec = r1_coarse_temp(:);
                A2_coarse =  finiteDiffGrid(2^(n-2)+1);
                e2_c2_jacobi = jacobi(A2_coarse, spdiags(diag(full(A2_coarse)), 0, size(A2_coarse,1),  size(A2_coarse,1)), zeros(size(r1_coarse_vec,1),1) , r1_coarse_vec, Njac);
%                 e2_c2_jacobi = A2_coarse\r1_coarse_vec;
                e2_c2_grid = reshape(e2_c2_jacobi, 2^(n-2)+1, 2^(n-2)+1)';
                
                 H = 1/(2^(n-2)+1);
                x=(0:2^(n-2))*H;
                y=(0:2^(n-2))*H;
                size(x)
                [x, y] = meshgrid(x,y);
        %         size(e_c_grid)
                surf(x, y,  e2_c2_grid , 'FaceColor','b')
                
                fprintf("Error and residual at grid level 3 ------- %d, %d \n", max(abs(e2_c2_jacobi )), max(abs(A2_coarse*e2_c2_jacobi-r1_coarse_vec)))  
                        % calculate residual - prepare to move to grid3
                        r2 = A2_coarse*e2_c2_jacobi - r1_coarse_vec;
                        r2_fine_grid = reshape(r2, 2^(n-2)+1, 2^(n-2)+1);
                        
                         % grid 4
                        r2_coarse = restriction(r2_fine_grid, n-3);
                        r2_coarse_temp = r2_coarse';
                        r2_coarse_vec = r2_coarse_temp(:);
                        A3_coarse =  finiteDiffGrid(2^(n-3)+1);
                        e3_c3_jacobi = jacobi(A3_coarse, spdiags(diag(full(A3_coarse)), 0, size(A3_coarse,1),  size(A3_coarse,1)), zeros(size(r2_coarse_vec,1),1) , r2_coarse_vec, Njac);
%                         e3_c3_jacobi = A3_coarse\r2_coarse_vec;
                        e3_c3_grid = reshape(e3_c3_jacobi, 2^(n-3)+1, 2^(n-3)+1)';
                         fprintf("Error and residual at grid level 4 ------- %d, %d \n", max(abs(e3_c3_jacobi )), max(abs(A3_coarse*e3_c3_jacobi-r2_coarse_vec)))  
                        
                e2_c2_grid_fine = interpolation(e3_c3_grid, n-2);
                e2_fine_temp = e2_c2_grid';
                e2_fine_vec = e2_fine_temp(:);
                e2_c2_jacobi_correction = e2_c2_jacobi + e2_fine_vec;
                e2_c2_jacobi_after_correction = jacobi(A2_coarse, spdiags(diag(full(A2_coarse)), 0, size(A2_coarse,1),  size(A2_coarse,1)),  e2_c2_jacobi_correction, r1_coarse_vec, Njac);
                e2_c2_grid_correction = reshape(e2_c2_jacobi_after_correction, 2^(n-2)+1, 2^(n-2)+1)' ;          
                        
                 fprintf("Error and residual at grid level 3 after correction ------- %d, %d\n ", max(abs(e2_c2_jacobi_after_correction )), max(abs(A2_coarse*e2_c2_jacobi_after_correction-r1_coarse_vec))) 
                
        e1_c1_grid = interpolation(e2_c2_grid_correction, n-1); 
%         e1_c1_grid = interpolation(e2_c2_grid, n-1);
        e1_fine_temp = e1_c1_grid';
        e1_fine_vec = e1_fine_temp(:);
        e_c_jacobi_correction = e_c_jacobi + e1_fine_vec;
        e_c_jacobi_after_correction = jacobi(A_coarse, spdiags(diag(full(A_coarse)), 0, size(A_coarse,1),  size(A_coarse,1)),  e_c_jacobi_correction, r_coarse_vec, Njac);
        e_c_grid_correction = reshape(e_c_jacobi_after_correction, 2^(n-1)+1, 2^(n-1)+1)';

        fprintf("Error and residual at grid level 2 after correction ------- %d, %d\n ", max(abs(e_c_jacobi_after_correction )), max(abs(A_coarse*e_c_jacobi_after_correction-r_coarse_vec)))  

    % step 5 : Interpolate the coarse grid error to the find grid
    % fine_grid_error = interpolation(e_c_grid, n);
    fine_grid_error = interpolation(e_c_grid_correction, n);
    % temp array becase I need the row major ordering in the vector
    e_fine_temp = fine_grid_error';
    e_fine_vec = e_fine_temp(:);
    
    % step 6 : correct the fine grid approximation
    ujacobi = ujacobi + e_fine_vec;
    
    % relax again with improved guess
    ujacobi = jacobi(A, D, ujacobi, f, Njac);
    u = ujacobi;
    fprintf("Error and residual at grid level 1 after correction ------- %d, %d \n", max(abs(ufd+ujacobi )), max(abs(A*ujacobi-f)))  
    
end

% from fine grid to coarse grid
function coarse = restriction(fine, n)
    % this only works if the fine grid has even number of rows and cols
    
    R = 2^n+1;
    C = R;
    
    coarse = zeros(R, C);
    coarse(1:end, 1:end) = fine(1:2:end, 1:2:end);
    
end

% interpolation from coarse grid to fine grid
% inputs - coarse grid, fine grid dimension
function fine = interpolation(v, n)

     R = 2^n+1;
    C = R;
    fine = zeros(R, C);
    
    fine(1:2:end, 1:2:end) = v;
    % rows
    fine(1:2:end, 2:2:end) = (fine(1:2:end, 3:2:end) + fine(1:2:end, 1:2:end-1)) /2;
    % columns
    fine(2:2:end, :) = (fine(3:2:end, :) + fine(1:2:end-1, :)) /2;

end

function Ak = finiteDiffGrid(N)
    Nt = N*N;
    h = 1/N;
    A = sparse(Nt,Nt);
    fact = 1/h/h;

    for i=1:N
        for j=1:N%N-1
            gi = (j-1)*N + i;
            A(gi,gi) = -4*fact;
        end
        for j=1:N-1%N-2
            gi = (j-1)*N + i;
            A(gi,gi+N) = fact;
        end
        for j=2:N%N-1
            gi = (j-1)*N + i;
            A(gi,gi-N) = fact;
        end
    end

    for j=1:N%N-1
        for i=1:(N-1)
            gi = (j-1)*N + i;
            A(gi,gi+1) = fact;
        end
        gi=j*N;
        A(gi,gi-N+1)= fact;
        for i=2:N
            gi = (j-1)*N + i;
            A(gi,gi-1) = fact;
        end
        gi=(j-1)*N+1;
        A(gi,gi+N-1)= fact;
    end
    Ak  = A;
end