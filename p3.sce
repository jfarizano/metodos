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

// Ej 2 Agrego el método de regula falsi para poder probarla
printf("\n\n---------Ej 2---------\n");
f1 = "cos(x) * cosh(x) + 1";
printf("Resultados f1(x) = %s\n", f1);
[y, i] = biseccion(f1, 1, 3, 10^(-6), 100);
printf("Bisección raiz %f en %d iteraciones\n", y, i);
[y, i] = secante(f1, 1, 3, 10^(-6), 100);
printf("Secante raiz %f en %d iteraciones\n", y, i);
[y, i] = regulaFalsi(f1, 1, 3, 10^(-6), 100);
printf("Regula falsi raiz %f en %d iteraciones\n", y, i);
printf("\n")

f2 = "2 * sin(x) - x^2"
printf("Resultados f2(x) = %s\n", f1);
[y, i] = biseccion(f2, 1, 2, 10^(-6), 100);
printf("Bisección raiz %f en %d iteraciones\n", y, i);
[y, i] = secante(f2, 1, 2, 10^(-6), 100);
printf("Secante raiz %f en %d iteraciones\n", y, i);
[y, i] = regulaFalsi(f2, 1, 2, 10^(-6), 100);
printf("Regula falsi raiz %f en %d iteraciones\n", y, i);
printf("\n")

f3 = "%e^(-x) - x^(-4)"
printf("Resultados f3(x) = %s\n", f1);
[y, i] = biseccion(f3, 1, 1.8, 10^(-6), 100);
printf("Bisección raiz %f en %d iteraciones\n", y, i);
[y, i] = secante(f3, 1, 1.8, 10^(-6), 100);
printf("Secante raiz %f en %d iteraciones\n", y, i);
[y, i] = regulaFalsi(f3, 1, 1.8, 10^(-6), 100);
printf("Regula falsi raiz %f en %d iteraciones\n", y, i);
printf("\n")

f4 = "log10(x) - x + 1"
printf("Resultados f4(x) = %s\n", f1);
[y, i] = biseccion(f4, 0.2, 1.4, 10^(-6), 100);
printf("Bisección raiz %f en %d iteraciones\n", y, i);
[y, i] = secante(f4, 0.2, 1.4, 10^(-6), 100);
printf("Secante raiz %f en %d iteraciones\n", y, i);
[y, i] = regulaFalsi(f4, 0.2, 1.4, 10^(-6), 100);
printf("Regula falsi raiz %f en %d iteraciones\n", y, i);
printf("\n")

f5 = "x^2 / 4 - sin(x)"
printf("Resultados f5(x) = %s\n", f1);
[y, i] = biseccion(f5, 1, 3, 10^(-6), 100);
printf("Bisección raiz %f en %d iteraciones\n", y, i);
[y, i] = secante(f5, 1, 3, 10^(-6), 100);
printf("Secante raiz %f en %d iteraciones\n", y, i);
[y, i] = regulaFalsi(f5, 1, 3, 10^(-6), 100);
printf("Regula falsi raiz %f en %d iteraciones\n", y, i);
printf("\n")

// Ej 5
printf("\n\n---------Ej 5---------\n");
g = "x + (1 / 5)*(x^2 - 5)";
[y, i] = puntoFijo(g, 2, 10^(-6), 100);
printf("g(x) = %s, y = %f, i = %d", g, y, i);

// Ej 6
printf("\n\n---------Ej 6---------\n");
g1 = "%e^x / 3"
[y, i] = puntoFijo(g1, 1, 10^(-6), 100);
printf("g1(x) = %s, y = %f, i = %d\n", g1, y, i)
g2 = "(%e^x - x) / 2"
[y, i] = puntoFijo(g2, 1, 10^(-6), 100);
printf("g2(x) = %s, y = %f, i = %d\n", g2, y, i)
g3 = "log(3 * x)"
[y, i] = puntoFijo(g3, 1, 10^(-6), 100);
printf("g3(x) = %s, y = %f, i = %d\n", g3, y, i)
g4 = "%e^x - 2 * x"
[y, i] = puntoFijo(g4, 1, 10^(-6), 100);
printf("g4(x) = %s, y = %f, i = %d\n", g4, y, i)


// Ej 7
printf("\n\n---------Ej 7---------\n");
function y = funEj7(x)
    y1 = 1 + x(1)^2 - x(2)^2 + exp(x(1)) * cos(x(2))
    y2 = 2 * x(1) * x(2) + exp(x(1)) * sin(x(2))
    y = [y1 ; y2]
endfunction

[yej7, iej7] = newtonMultivariable(funEj7, [-1 4], 10^(-12), 5)

printf("y = [%f %f], i = %d", yej7(1), yej7(2), iej7)

// Ej 8
printf("\n\n---------Ej 8---------\n");
function y = funEj8(x)
    y1 = x(1)^2 + x(1) * (x(2)^3) - 9
    y2 = 3 * (x(1)^2) * x(2) - 4 - x(2)^3
    y = [y1 ; y2]
endfunction

xej8 = [1.2 2.5]
[yej8, iej8] = newtonMultivariable(funEj8, xej8, 10^(-6), 100)
printf("x = [%f %f], y = [%f %f], i = %d\n", xej8(1), xej8(2), yej8(1), yej8(2), iej8)

xej8 = [-2 2.5]
[yej8, iej8] = newtonMultivariable(funEj8, xej8, 10^(-6), 100)
printf("x = [%f %f], y = [%f %f], i = %d\n", xej8(1), xej8(2), yej8(1), yej8(2), iej8)

xej8 = [-1.2 -2.5]
[yej8, iej8] = newtonMultivariable(funEj8, xej8, 10^(-6), 100)
printf("x = [%f %f], y = [%f %f], i = %d\n", xej8(1), xej8(2), yej8(1), yej8(2), iej8)

xej8 = [2 -2.5]
[yej8, iej8] = newtonMultivariable(funEj8, xej8, 10^(-6), 100)
printf("x = [%f %f], y = [%f %f], i = %d\n", xej8(1), xej8(2), yej8(1), yej8(2), iej8)

// Ej 9

// Ej 10
printf("\n\n---------Ej 10---------\n");
printf("Punto fijo empezando desde x0 = 28\n")
// Acá busco una f(x) tal que f(x) = x
f = "9.8 * 25 / (2 * %pi) * tanh(8 * %pi / x)"
[y, i] = puntoFijo(f, 28, 10^(-1), 1000);
printf("f(x) = %s\n y = %f, i = %d\n", f, y, i)

printf("Con el método de newton, con x0 = %f obtenido del cálculo anterior\n", y)
// Busco g(x) tal que g(x) = 0, entonces g(x) = f(x) - x
[y, i] = newton(f + " - x", y, 10^(-6), 10000, 10^(-6))
printf("f(x) = %s\n y = %f, i = %d\n", f, y, i)
