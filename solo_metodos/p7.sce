exec("solo_metodos/p4.sce", -1) // Necesito eliminación gaussiana

// Diferencias dividas
function w=DD(x,y)
    n = length(x);
    if n==1 then
        w = y(1)
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end;
endfunction

// Polinomio interpolante (con Newton)
function p = interpolNewton(x,y)
    r = poly(0,"x");
    p = 0;
    n= length(x);
    for i=n:(-1):2
        p = (p+DD(x(1:i),y(1:i)))*(r-x(i-1))
    end;
    p = p + y(1);
endfunction

// Interpolación Lagrange
function y = Lk(x,k)
    n = length(x)
    r = [x(1:k-1) x(k+1:n)]
    p = poly(r,"x","roots")
    pk = horner(p,x(k))
    y = p / pk
endfunction

function z = interpolLagrange(x,y)
    n = length(x)
    pol = 0
    for k = 1:n
        pol = pol + (Lk(x,k)*y(k))
    end
    z = pol
endfunction

// Aproximación polinomial de mínimos cuadrados para matrices con rango completo
function [p,err] = minCuadPol(A,b)
    [a,A_amp] = gaussElimPP((A')*A,(A')*(b'))
    p = poly(a,"x","coeff")
    err = norm(A*a-(b'))
endfunction

function [p,err] = minCuadPolInv(A,b)
    a = inv(A'* A) * A' * b'
    p = poly(a,"x","coeff")
    err = norm(A*a-(b'))
endfunction

function [p,err] = minCuadPolQR(A,b)
    [Q, R] = qr(A, ["e"])
    a = inv(R) * Q' * b'
    p = poly(a,"x","coeff")
    err = norm(A*a-(b'))
endfunction

// Matriz del método de mínimo cuadrados polinomial
function A = matrizMinCuadPol(x,grado)
    // p = grado+1
    m = length(x)
    A = ones(m,1)
    for j=1:grado
        A = [A,(x').^j]
    end
endfunction

// Raíces del polinomio de Chebyshev
function w = rootsCheby(n)
    // Entrada: n = grad del polinomio de Chebyshev
    // Salida: w = vector con las raices del polinomio de Chebyshev
    for i=0:(n-1) do
        w(i+1)=cos((2*i+1)*%pi/(2*n))
    end
    w = w'
endfunction

// Raíces del polinomio de Chebyshev
function w = rootsCheby_ab(n,a,b)
    for i=0:(n-1) do
        w(i+1)=cos((2*i+1)*%pi/(2*n))
    end
    w = ((b+a) + w*(b-a))/2
    w = w'
endfunction

// Polinomio de Chebyshev
function w = polCheby(x,n)
    // Entrada: n = número natural; x = número real
    // Salida: Polinomio de Chebyshev de grado n evaluado en x
    if n==0 then
        w = 1
    elseif n==1 then
        w = x
    elseif n==2 then
        w = 2*x.^2-1
    else
        w = 2*x.*polCheby(x,n-1)-polCheby(x,n-2)
    end
endfunction

// Funcion utilizada para acotar más rápidamente el error de una interpolación
// tener en cuenta que estoy multiplicando por el máximo teórico de la derivada
// calculado en papel, pero sirve para graficar el polinomio y ver el máximo
// para encontrar la cota
function p = polErrorInterpol(x, max_deriv)
    n = size(x, 2)

    x_pol = poly(0, "x")
    p = 1

    for i=1:n
        p = p * (x_pol - x(i))
    end
    p = p * (max_deriv / factorial(n))
endfunction
