chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/errores.sce", -1)
exec("solo_metodos/p7.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
x = [-2.0 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2.0]
y = [1.50 0.99 0.61 0.27 0.02 -0.0096 -0.065 0.38 0.63 0.98 1.50]
printf("x = ")
disp(x)
printf("y = ")
disp(y)
printf("Matriz del método")
A = matrizMinCuadPol(x, 4)
disp(A)
printf("Polinomio")
[p, err] = minCuadPol(A, y)
disp(p)
xx = -2:0.001:2
plot2d(x', y, 1, leg="tabla")
plot2d(xx', horner(p, xx)', 2, leg="p4(x)")


// Ejercicio 2
printf("\n\n---------Ej 2---------\n");
clf()
printf("f(x) = e^x\n")
printf("Nodos de interpolación equidistantes\n")
x = [-1 -1/3 1/3 1]
y = exp(x)
printf("x = ")
disp(x)
printf("y = ")
disp(y)
p = interpolNewton(x, y)
printf("Polinomio interpolante con nodos equidistantes")
disp(p)
printf("\ne^(-0.09)\n")
imprimir_errores(exp(-0.09), horner(p, -0.09))
printf("\ne^(0.01)\n")
imprimir_errores(exp(0.01), horner(p, 0.01))


printf("\nNodos de Chebyshev\n")
x_cheby = rootsCheby(4)
y_cheby = exp(x_cheby)
printf("x = ")
disp(x_cheby)
printf("y = ")
disp(y_cheby)
q = interpolNewton(x_cheby, y_cheby)
printf("Polinomio interpolante con nodos de Chebyshev")
disp(q)
printf("\ne^(-0.09)\n")
imprimir_errores(exp(-0.09), horner(q, -0.09))
printf("\ne^(0.01)\n")
imprimir_errores(exp(0.01), horner(q, 0.01))

xx = -1:0.001:1
e_p = exp(xx) - horner(p, xx)
e_q = exp(xx) - horner(q, xx)
plot2d(xx', [e_p' e_q'], [2, 3], leg = "p(x)@q(x)")

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
clf()
x = 1:10
y = [32.9 30.8 26.4 24.2 19.2 16.5 19.3 21 23 26.2]
printf("x = ")
disp(x)
printf("y = ")
disp(y)
grados = [3 5 7 9]
xx = 1:0.1:10
yy = []
leyen = "tabla"

for i = 1:size(grados, 2)
    grado = grados(i)
    A = matrizMinCuadPol(x, grado)
    [p, err] = minCuadPolInv(A, y)
    yy = [yy horner(p, xx')]
    leyen = leyen + ("@p"+string(grado)+"(x)")
end

plot2d(x, y)
plot2d(xx', yy, 1:size(grados, 2), leg=leyen)

clf()
yy = []
for i = 1:size(grados, 2)
    grado = grados(i)
    A = matrizMinCuadPol(x, grado)
    [p, err] = minCuadPolQR(A, y)
    yy = [yy horner(p, xx')]
    leyen = leyen + ("@p"+string(grado)+"(x)")
end

plot2d(x, y)
plot2d(xx', yy, 1:size(grados, 2), leg=leyen)

// Salida completa del programa
/*
---------Ej 1---------
x = 
  -2.  -1.6  -1.2  -0.8  -0.4   0.   0.4   0.8   1.2   1.6   2.
y = 
   1.5   0.99   0.61   0.27   0.02  -0.0096  -0.065   0.38   0.63   0.98   1.5
Matriz del método
   1.  -2.    4.    -8.      16.   
   1.  -1.6   2.56  -4.096   6.5536
   1.  -1.2   1.44  -1.728   2.0736
   1.  -0.8   0.64  -0.512   0.4096
   1.  -0.4   0.16  -0.064   0.0256
   1.   0.    0.     0.      0.    
   1.   0.4   0.16   0.064   0.0256
   1.   0.8   0.64   0.512   0.4096
   1.   1.2   1.44   1.728   2.0736
   1.   1.6   2.56   4.096   6.5536
   1.   2.    4.     8.      16.   
Polinomio
                                                     2                  3                 4
  -0.034551981352 +0.0183906371406x +0.4778809731935x  -0.0052204739705x  -0.024443655303x 


---------Ej 2---------
f(x) = e^x
Nodos de interpolación equidistantes
x = 
  -1.  -0.3333333333333   0.3333333333333   1.
y = 
   0.3678794411714   0.7165313105738   1.3956124250861   2.718281828459
Polinomio interpolante con nodos equidistantes
                                                     2                  3
   0.9951957719568 +0.999049231534x +0.5478848628585x  +0.1761519621098x 

e^(-0.09)
Valor real = 0.913931 ~ Valor aproximado = 0.909591
Error absoluto =
   0.0043403915437
Error relativo =
   0.0047491448084

e^(0.01)
Valor real = 1.010050 ~ Valor aproximado = 1.005241
Error absoluto =
   0.0048089381738
Error relativo =
   0.0047610884395

Nodos de Chebyshev
x = 
   0.9238795325113   0.3826834323651  -0.3826834323651  -0.9238795325113
y = 
   2.519044171407   1.4662138007571   0.6820287733505   0.3969759686435
Polinomio interpolante con nodos de Chebyshev
                                                     2                  3
   0.994615316879 +0.9989332279763x +0.5429007233211x  +0.1751756940472x 

e^(-0.09)
Valor real = 0.913931 ~ Valor aproximado = 0.908981
Error absoluto =
   0.0049500661322
Error relativo =
   0.0054162350645

e^(0.01)
Valor real = 1.010050 ~ Valor aproximado = 1.004659
Error absoluto =
   0.0053910526774
Error relativo =
   0.005337410807


---------Ej 4---------
x = 
   1.   2.   3.   4.   5.   6.   7.   8.   9.   10.
y = 
   32.9   30.8   26.4   24.2   19.2   16.5   19.3   21.   23.   26.2
*/
