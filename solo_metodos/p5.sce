function [x1, contador] = jacobi(A, b, x0, tol, iter)
    n = size(A, 1)
    x1 = x0

    contador = 0
    delta = %inf

    while delta > tol && contador < iter
        x1(1) = (1/A(1,1)) * (b(1) - A(1,2:n)*x0(2:n))
        for i=2:n-1
            x1(i) = (1/A(i,i)) * (b(i) - A(i,1:i-1)*x0(1:i-1) - A(i,i+1:n)*x0(i+1:n))
        end
        x1(n) = (1/A(n,n)) * (b(n) - A(n,1:n-1)*x0(1:n-1))

        delta = norm(x1 - x0)
        x0 = x1
        contador = contador + 1
    end
endfunction

function [x1, contador] = gaussSeidel(A, b, x0, tol, iter)
    n = size(A, 1)
    x1 = x0

    contador = 0
    delta = %inf

    while delta > tol && contador < iter
        x1(1) = (1/A(1,1)) * (b(1) - A(1,2:n)*x1(2:n))
        for i=2:n-1
            x1(i) = (1/A(i,i)) * (b(i) - A(i,1:i-1)*x1(1:i-1) - A(i,i+1:n)*x1(i+1:n))
        end
        x1(n) = (1/A(n,n)) * (b(n) - A(n,1:n-1)*x1(1:n-1))

        delta = norm(x1 - x0)
        x0 = x1
        contador = contador + 1
    end
endfunction

function w = omega_opt(A)
    // A es una matriz definida positiva y tridiagonal
    I = eye(A)
    invD = diag(1 ./ diag(A))
    Tj = (I - invD * A)
    // El máximo de los absolutos de los autovalores
    rho = max(abs(spec(Tj)))

    w = 2/(1+sqrt(1 - rho**2))
endfunction

function [x1, contador] = SOR(A, b, x0, w, tol, iter)
    n = size(A, 1)
    x1 = x0

    contador = 0
    delta = %inf

    while delta > tol && contador < iter
        x1(1) = (1-w)*x0(1) + (w/A(1,1)) * (b(1) - A(1,2:n)*x1(2:n))
        for i=2:n-1
            x1(i) = (1-w)*x0(i) + (w/A(i,i)) * (b(i) - A(i,1:i-1)*x1(1:i-1) - A(i,i+1:n)*x1(i+1:n))
        end
        x1(n) = (1-w)*x0(n) + (w/A(n,n)) * (b(n) - A(n,1:n-1)*x1(1:n-1))

        delta = norm(x1 - x0)
        x0 = x1
        contador = contador + 1
    end
endfunction

// Funciones para obtener el radio espectral de las matrices de los métodos
// de Jacobi, Gauss-Seidel y SOR para verificar su convergencia
function r = r_Jacobi(A)
   D_inv = diag(diag(1 ./ A))

   T = eye(A) - (D_inv)
   r = max(abs(spec(T)))
endfunction

function r = r_GaussSeidel(A)
   n = size(A, 1)
   D = diag(diag(A))
   L = A - D

   for i=1:n
       for j=1:n
           if i > j
               L(i, j) = 0
           end
       end
   end

   N = L + D
   T = eye(A) - inv(N)*A
   r = max(abs(spec(T)))
endfunction

// El radio depende del omega utilizado
function r = r_SOR(A, w)
   n = size(A, 1)
   D = diag(diag(A))
   L = A - D
   U = A - D

   for i=1:n
       for j=1:n
           if i > j
               L(i, j) = 0
           end
           if i < j
               U(i, j) = 0
           end
       end
   end

   Tw = inv((D + w*L)) * (((1-w) * D) - (w*U))
   r = max(abs(spec(Tw)))
endfunction
