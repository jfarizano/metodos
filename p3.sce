chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p3.sce", -1)

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

// Salida completa del programa
/*
---------Ej 2---------
Resultados f1(x) = cos(x) * cosh(x) + 1
Bisección raiz 1.875104 en 21 iteraciones
Secante raiz 1.875104 en 9 iteraciones
Regula falsi raiz 1.875104 en 22 iteraciones

Resultados f2(x) = cos(x) * cosh(x) + 1
Bisección raiz 1.404415 en 18 iteraciones
Secante raiz 1.404415 en 7 iteraciones
Regula falsi raiz 1.404415 en 13 iteraciones

Resultados f3(x) = cos(x) * cosh(x) + 1
Bisección raiz 1.429611 en 18 iteraciones
Secante raiz 1.429612 en 13 iteraciones
Regula falsi raiz 1.429616 en 34 iteraciones

Resultados f4(x) = cos(x) * cosh(x) + 1
Bisección raiz 1.000002 en 18 iteraciones
Secante raiz 1.000000 en 7 iteraciones
Regula falsi raiz 0.999999 en 7 iteraciones

Resultados f5(x) = cos(x) * cosh(x) + 1
Bisección raiz 1.933754 en 18 iteraciones
Secante raiz 1.933754 en 8 iteraciones
Regula falsi raiz 1.933753 en 14 iteraciones



---------Ej 5---------
g(x) = x + (1 / 5)*(x^2 - 5), y = -2.236068, i = 14

---------Ej 6---------
g1(x) = %e^x / 3, y = 0.619062, i = 28
g2(x) = (%e^x - x) / 2, y = 0.619062, i = 17
g3(x) = log(3 * x), y = 1.512133, i = 32
g4(x) = %e^x - 2 * x, y = 0.619061, i = 8


---------Ej 7---------

 Se alcanzó el máximo de iteraciones
y = [-0.293163 1.172660], i = 5

---------Ej 8---------
x = [1.200000 2.500000], y = [1.336355 1.754235], i = 5
x = [-2.000000 2.500000], y = [-0.901266 -2.086588], i = 9
x = [-1.200000 -2.500000], y = [-0.901266 -2.086588], i = 5
x = [2.000000 -2.500000], y = [-3.001625 0.148108], i = 19


---------Ej 10---------
Punto fijo empezando desde x0 = 28
f(x) = 9.8 * 25 / (2 * %pi) * tanh(8 * %pi / x)
 y = 27.955296, i = 2
Con el método de newton, con x0 = 27.955296 obtenido del cálculo anterior
f(x) = 9.8 * 25 / (2 * %pi) * tanh(8 * %pi / x)
 y = 27.928560, i = 3
*/
