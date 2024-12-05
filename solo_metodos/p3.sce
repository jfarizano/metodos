function [c, i] = biseccion(fun, a, b, tol, iter)
    deff("y=f(x)", "y="+fun);
    i = 0;
    delta = %inf;

    while delta > tol && i < iter
        i = i + 1;

        c = (a + b) / 2;
        fc = f(c);

        if abs(fc) < tol then
            return;
        end

        delta = b - c;
        if f(a) * fc < 0 then
            b = c;
        else
            a = c;
        end
    end

    if (delta > tol) then c = %nan; disp('Se alcanzó el máximo de iteraciones')
    end
endfunction

function [x1, i] = secante(fun, x0, x1, tol, iter)
    deff("y=f(x)", "y="+fun);
    i = 0;

    // Utilizo x0 para la iteración anterior y x1 la actual
    fx0 = f(x0)
    fx1 = f(x1)
    delta = abs(fx1 - fx0)

    while delta > tol && i < iter
        // Calculo la próxima iteración
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        // Ciclo las variables para que en el proximo loop la nueva iteración sea la
        // actual y la actual pase a ser la anterior
        x0 = x1
        fx0 = fx1
        x1 = x2
        fx1 = f(x1)

        delta = abs(fx1 - fx0)
        i = i + 1
    end

    if (delta > tol) then x1 = %nan; disp('Se alcanzó el máximo de iteraciones')
    end
endfunction

function [x1, i] = newton(fun, x0, tol, iter, h)
    deff("y=f(x)", "y="+fun);
    i = 0;
    delta = %inf;

    while delta > tol && i < iter
        fx0 = f(x0);
        dfx0 = f(x0+h) - fx0;
        x1 = x0 - fx0 * h / dfx0;

        delta = abs(x1 - x0);
        x0 = x1;
        i = i + 1;
    end

    if (delta > tol) then x1 = %nan; disp('Se alcanzó el máximo de iteraciones')
    end
endfunction

function [c, i] = regulaFalsi(fun, a, b, tol, iter)
    deff("y=f(x)", "y="+fun);
    i = 0;
    delta = %inf;

    while delta > tol && i < iter
        i = i + 1;

        fa = f(a);
        fb = f(b);
        c = b - fb * ((b - a) / (fb - fa));
        fc = f(c);

        if abs(fc) < tol then
            return;
        end

        delta = abs(fc - fb);
        if fa * fc < 0 then
            b = c;
        else
            a = c;
        end
    end

    if (delta > tol) then c = %nan; disp('Se alcanzó el máximo de iteraciones')
    end
endfunction

function [x1, i] = puntoFijo(fun, x0, tol, iter)
    deff("y=f(x)", "y="+fun);
    i = 0;
    delta = %inf;

    while delta > tol && i < iter
        x1 = f(x0);

        delta = abs(x1 - x0);
        x0 = x1;
        i = i + 1;
    end

    if (delta > tol) then x1 = %nan; disp('Se alcanzó el máximo de iteraciones')
    end
endfunction

// fun es una función de Scilab que recibe un vector x, accede a los x(i)
// y devuelve un vector "y" con la solución
function [x1, i] = newtonMultivariable(fun, x0, tol, iter)
    i = 0;
    delta = %inf;
    x0 = x0';

    while delta > tol && i < iter
        J = numderivative(fun, x0);
        invJ = inv(J);
        x1 = x0 - (invJ * fun(x0));

        delta = norm(x1 - x0);
        x0 = x1;
        i = i + 1;
    end

    if (delta > tol) then disp('Se alcanzó el máximo de iteraciones')
    end
endfunction
