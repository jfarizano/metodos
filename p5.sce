chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p4.sce", -1) // Importo p4 para traer eliminación de gaus con pp
exec("solo_metodos/p5.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
A = [0  2  4;
     1 -1 -1;
     1 -1  2]
b = [0 0.375 0]'
x0 = [0 0 0]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x0 = \n")
disp(x0)
//x1 = jacobi(A, b, x0, 10^(-6), %inf)
//printf("Con método de Jacobi: x = \n")
//disp(x1)
//x1 = gaussSeidel(A, b, x0, 10^(-6), %inf)
//printf("Con método de Gauss-Seidel: x = \n")
//disp(x1)
printf("No convergen\n")

A = [ 1 -1  0;
     -1  2 -1;
      0 -1  1.1]
b = [0 1 0]'
x0 = [0 0 0]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x0 = \n")
disp(x0)
x1 = jacobi(A, b, x0, 10^(-6), %inf)
printf("Con método de Jacobi: x = \n")
disp(x1)
x1 = gaussSeidel(A, b, x0, 10^(-6), %inf)
printf("Con método de Gauss-Seidel: x = \n")
disp(x1)

// Ejercicio 2
printf("\n\n---------Ej 2---------\n");
A = [10 1  2  3   4;
     1  9 -1  2  -3;
     2 -1  7  3  -5;
     3  2  3  12 -1;
     4 -3 -5 -1   15]
b = [12 -27 14 -17 12]'
x0 = [0 0 0 0 0]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x0 = \n")
disp(x0)
[x1, iters] = jacobi(A, b, x0, 10^(-6), %inf)
printf("Con método de Jacobi: x = \n")
disp(x1)
printf("En %d iteraciones\n", iters)
[x1, iters] = gaussSeidel(A, b, x0, 10^(-6), %inf)
printf("Con método de Gauss-Seidel: x = \n")
disp(x1)
printf("En %d iteraciones\n", iters)

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
function [A, b] = matrizEj4(N)
    A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3)
    b = ones(N,1)
endfunction

N = 100
[A, b] = matrizEj4(N)
printf("N = %d\n", N)
tic()
x = gaussElimPP(A, b)
t = toc()
printf("Tiempo con eliminación de Gauss con pivoteo = %f segundos\n",t)
printf("Método de Gauss-Seidel:\n")
x0 = zeros(N, 1)
tic()
[x1, iters] = gaussSeidel(A, b, x0, 10^(-6), %inf)
t = toc()
printf("Tolerancia de 10^(-6): %f segundos y %d iteraciones\n", t, iters)
tic()
[x1, iters] = gaussSeidel(A, b, x0, 10^(-11), %inf)
t = toc()
printf("Tolerancia de 10^(-11): %f segundos y %d iteraciones\n", t, iters)
printf("\n\n")

N = 500
[A, b] = matrizEj4(N)
printf("N = %d\n", N)
tic()
x = gaussElimPP(A, b)
t = toc()
printf("Tiempo con eliminación de Gauss con pivoteo = %f segundos\n",t)
printf("Método de Gauss-Seidel:\n")
x0 = zeros(N, 1)
tic()
[x1, iters] = gaussSeidel(A, b, x0, 10^(-6), %inf)
t = toc()
printf("Tolerancia de 10^(-6): %f segundos y %d iteraciones\n", t, iters)
[x1, iters] = gaussSeidel(A, b, x0, 10^(-11), %inf)
t = toc()
printf("Tolerancia de 10^(-11): %f segundos y %d iteraciones\n", t, iters)
printf("\n\n")


// Ejercicio 5
printf("\n\n---------Ej 5---------\n");
A = [4  3  0;
     3  4 -1;
     0 -1  4]
b = [24 30 -24]'
x0 = [0 0 0]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x0 = \n")
w = 0.5 // omega arbitrario
[x1, iters] = SOR(A, b, x0, w, 10^(-2), %inf)
printf("Aprox inicial con SOR con omega = %f\n", w)
disp(x1)
printf("En %d iteraciones\n", iters)
x0 = x1 // Pasa a ser el vector inicial de las siguientes aproximaciones

// Item a)
[x1, iters] = gaussSeidel(A, b, x0, 10^(-7), %inf)
printf("Aproximación con método de Gauss-Seidel, x = \n")
disp(x1)
printf("En %d iteraciones\n", iters)

// Item b)
w = omega_opt(A)
printf("Con método de SOR y omega óptimo = %f\n", w)
[x1, iters] = SOR(A, b, x0, w, 10^(-7), %inf)
printf("x = \n")
disp(x1)
printf("En %d iteraciones\n", iters)

// Salida completa del programa
/*
---------Ej 1---------
A = 

   0.   2.   4.
   1.  -1.  -1.
   1.  -1.   2.
b = 

   0.
   0.375
   0.
x0 = 

   0.
   0.
   0.
No convergen
A = 

   1.  -1.   0. 
  -1.   2.  -1. 
   0.  -1.   1.1
b = 

   0.
   1.
   0.
x0 = 

   0.
   0.
   0.
Con método de Jacobi: x = 

   10.999978921756
   10.999979879858
   9.9999808379602
Con método de Gauss-Seidel: x = 

   10.999987364362
   10.999987938709
   9.9999890351902


---------Ej 2---------
A = 

   10.   1.   2.   3.    4. 
   1.    9.  -1.   2.   -3. 
   2.   -1.   7.   3.   -5. 
   3.    2.   3.   12.  -1. 
   4.   -3.  -5.  -1.    15.
b = 

   12.
  -27.
   14.
  -17.
   12.
x0 = 

   0.
   0.
   0.
   0.
   0.
Con método de Jacobi: x = 

   1.0000016240578
  -2.0000014882294
   2.9999972712683
  -1.9999995623834
   0.999998051336
En 67 iteraciones
Con método de Gauss-Seidel: x = 

   1.0000009087571
  -2.0000007460373
   2.999998708944
  -1.9999998791567
   0.9999991861615
En 38 iteraciones


---------Ej 4---------
N = 100
Tiempo con eliminación de Gauss con pivoteo = 0.044316 segundos
Método de Gauss-Seidel:
Tolerancia de 10^(-6): 0.011058 segundos y 18 iteraciones
Tolerancia de 10^(-11): 0.022564 segundos y 40 iteraciones


N = 500
Tiempo con eliminación de Gauss con pivoteo = 1.936798 segundos
Método de Gauss-Seidel:
Tolerancia de 10^(-6): 0.096383 segundos y 18 iteraciones
Tolerancia de 10^(-11): 0.307122 segundos y 40 iteraciones




---------Ej 5---------
A = 

   4.   3.   0.
   3.   4.  -1.
   0.  -1.   4.
b = 

   24.
   30.
  -24.
x0 = 
Aprox inicial con SOR con omega = 0.500000

   3.040986903533
   3.9597852562084
  -5.0117816996175
En 16 iteraciones
Aproximación con método de Gauss-Seidel, x = 

   3.0000000945384
   3.999999921218
  -5.0000000196955
En 28 iteraciones
Con método de SOR y omega óptimo = 1.240408
x = 

   3.0000000096426
   3.9999999948158
  -5.0000000008117
En 13 iteraciones
*/
