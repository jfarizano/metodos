function [x1, contador] = jacobi(A, b, x0, tol, iter)
    n = size(A, 1)
    x1 = x0

    contador = 0
    delta = %inf

    while delta > tol && contador < iter
        for i=1:n
            suma = 0
            for j=1:n
                if (j<>i)
                    suma = suma + A(i,j) * x0(j)
                end
            end
            x1(i) = 1/A(i,i) * (b(i) - suma)
        end

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
        for i=1:n
            suma = 0
            for j=1:n
                if (j<>i)
                    suma = suma + A(i,j) * x1(j)
                end
            end
            x1(i) = 1/A(i,i) * (b(i) - suma)
        end

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
        for i=1:n
            suma = 0
            for j=1:n
                if (j<>i)
                    suma = suma + A(i,j) * x1(j)
                end
            end
            x1(i) = (1-w) * x1(i) + w * 1/A(i,i) * (b(i) - suma)
        end

        delta = norm(x1 - x0)
        x0 = x1
        contador = contador + 1
    end
endfunction
