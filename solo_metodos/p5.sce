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
