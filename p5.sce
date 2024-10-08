function x1 = jacobi(A, b, x0, eps)
    n = size(A, 1)
    x1 = x0

    cont = 0
    delta = %inf

    while(delta > eps)
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
        cont = cont + 1
    end

    disp(cont)
endfunction

function x1 = gaussSeidel(A, b, x0, eps)
    n = size(A, 1)
    x1 = x0

    cont = 0
    delta = %inf

    while(delta > eps)
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
        cont = cont + 1
    end

    disp(cont)
endfunction

A = [ 1 -1  0;
     -1  2 -1;
      0 -1  1.1]

b = [0 1 0]'
x0 = [0 0 0]'
x1 = jacobi(A, b, x0, 10^(-6))
disp(x1)
x1 = gaussSeidel(A, b, x0, 10^(-6))
disp(x1)
