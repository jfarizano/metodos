chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p6.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
A = [ 1  0 0;
     -1  0 1;
     -1 -1 2;]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
printf("\n\n")


A = [1    0   0;
    -0.1  0   0.1;
    -0.1 -0.1 2]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
printf("\n\n")

A = [1     0    0;
    -0.25  0    0.25
    -0.25 -0.25 2]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
printf("\n\n")

A = [4 -1  0;
    -1  4 -1;
    -1 -1  4]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
printf("\n\n")

A = [3 2 1;
     2 3 0;
     1 0 3]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
printf("\n\n")

A = [ 4.75 2.25 -0.25;
      2.25 4.74  1.25;
     -0.25 1.25  4.75]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))

// Ejercicio 3
printf("\n\n---------Ej 3---------\n");
function A = matrizEj3(eps)
A = [ 1  -1  0;
     -2   4 -2;
      0  -1  1+eps]
endfunction

for k=0:10
    eps = 0.1 * k
    printf("eps = %f\n", eps)
    A = matrizEj3(eps)
    printf("A = \n")
    disp(A)
    printf("Polinomio característico = \n")
    pol = poly(A, "x")
    disp(pol)
    printf("Raices del pol = \n")
    r = roots(pol)
    disp(r)
    printf("Autovalores = \n")
    disp(spec(A))
    printf("\n\n")
end

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
// Matriz de prueba sacada de la teoría
A = [1 -1;
     2 -1]
printf("A = \n")
disp(A)
circ = circulosGers(A)
printf("Cotas de autovalores según el teorema de Gerschorin\n")
for i=1:size(A, 1)
    centro = circ(i, 1)
    radio = circ(i, 2)
    a = centro - radio
    b = centro + radio
    printf("[%f, %f] ", a, b)
end
printf("\n")
printf("Autovalores según scilab: \n")
disp(spec(A))
dibujarCirculosGersValor(A)

// Ejercicio 5
printf("\n\n---------Ej 5---------\n");
A = [6 4 4 1;
     4 6 1 4;
     4 1 6 4;
     1 4 4 6]
printf("A = \n")
disp(A)
[valor, z, iters] = potencia(A, [1 1 0 0]', 10^(-8), %inf)
printf("Autovalor dominante = %f\n", valor)
printf("Autovector asociado = \n")
disp(z)
printf("Calculado en %d iteraciones\n", iters)
valor = maxAutovalor(A)
printf("Autovalor dominante según scilab = %f\n", valor)
printf("\n\n")

A = [12 1  3  4; 
     1 -3  1  5; 
     3  1  6 -2; 
     4  5 -2 -1]
printf("A = \n")
disp(A)
[valor, z, iters] = potencia(A, [1 1 0 0]', 10^(-8), %inf)
printf("Autovalor dominante = %f\n", valor)
printf("Autovector asociado = \n")
disp(z)
printf("Calculado en %d iteraciones\n", iters)
valor = maxAutovalor(A)
printf("Autovalor dominante según scilab = %f\n", valor)


// Salida completa del programa
/*
---------Ej 1---------
A = 

   1.   0.   0.
  -1.   0.   1.
  -1.  -1.   2.
Cotas de autovalores según el teorema de Gerschorin
[1.000000, 1.000000] [-2.000000, 2.000000] [0.000000, 4.000000] 
Autovalores según scilab: 

   1.  
   1.  
   1.  


A = 

   1.    0.    0. 
  -0.1   0.    0.1
  -0.1  -0.1   2. 
Cotas de autovalores según el teorema de Gerschorin
[1.000000, 1.000000] [-0.200000, 0.200000] [1.800000, 2.200000] 
Autovalores según scilab: 

   1.9949874371066  
   0.0050125628934  
   1.  


A = 

   1.     0.     0.  
  -0.25   0.     0.25
  -0.25  -0.25   2.  
Cotas de autovalores según el teorema de Gerschorin
[1.000000, 1.000000] [-0.500000, 0.500000] [1.500000, 2.500000] 
Autovalores según scilab: 

   1.9682458365519  
   0.0317541634481  
   1.  


A = 

   4.  -1.   0.
  -1.   4.  -1.
  -1.  -1.   4.
Cotas de autovalores según el teorema de Gerschorin
[3.000000, 5.000000] [2.000000, 6.000000] [2.000000, 6.000000] 
Autovalores según scilab: 

   4.6180339887499  
   2.3819660112501  
   5.  


A = 

   3.   2.   1.
   2.   3.   0.
   1.   0.   3.
Cotas de autovalores según el teorema de Gerschorin
[0.000000, 6.000000] [1.000000, 5.000000] [2.000000, 4.000000] 
Autovalores según scilab: 

   0.7639320225002
   3.
   5.2360679774998


A = 

   4.75   2.25  -0.25
   2.25   4.74   1.25
  -0.25   1.25   4.75
Cotas de autovalores según el teorema de Gerschorin
[2.250000, 7.250000] [1.240000, 8.240000] [3.250000, 6.250000] 
Autovalores según scilab: 

   2.0598495478351
   4.9616720796854
   7.2184783724795


---------Ej 3---------
eps = 0.000000
A = 

   1.  -1.   0.
  -2.   4.  -2.
   0.  -1.   1.
Polinomio característico = 

                          2   3
  -2.257727500D-16 +5x -6x  +x 
Raices del pol = 

   5.  
   1.  
   0.  
Autovalores = 

   5.  
   1.  
   0.  


eps = 0.100000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.1
Polinomio característico = 

                  2   3
  -0.2 +5.5x -6.1x  +x 
Raices del pol = 

   5.0102087889049  
   1.0518401043293  
   0.0379511067658  
Autovalores = 

   5.0102087889049  
   0.0379511067658  
   1.0518401043293  


eps = 0.200000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.2
Polinomio característico = 

                2   3
  -0.4 +6x -6.2x  +x 
Raices del pol = 

   5.0208507939808  
   1.1071945792249  
   0.0719546267943  
Autovalores = 

   5.0208507939808  
   0.0719546267943  
   1.1071945792249  


eps = 0.300000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.3
Polinomio característico = 

                  2   3
  -0.6 +6.5x -6.3x  +x 
Raices del pol = 

   5.0319505983321  
   1.1657664301585  
   0.1022829715094  
Autovalores = 

   5.0319505983321  
   0.1022829715094  
   1.1657664301585  


eps = 0.400000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.4
Polinomio característico = 

                2   3
  -0.8 +7x -6.4x  +x 
Raices del pol = 

   5.0435343569934  
   1.2272144557583  
   0.1292511872483  
Autovalores = 

   5.0435343569934  
   0.1292511872483  
   1.2272144557583  


eps = 0.500000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.5
Polinomio característico = 

                2   3
  -1 +7.5x -6.5x  +x 
Raices del pol = 

   5.0556298872046  
   1.2911771179572  
   0.1531929948382  
Autovalores = 

   5.0556298872046  
   0.1531929948382  
   1.2911771179572  


eps = 0.600000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.6
Polinomio característico = 

                2   3
  -1.2 +8x -6.6x  +x 
Raices del pol = 

   5.0682667593719  
   1.3572923111759  
   0.1744409294521  
Autovalores = 

   5.0682667593719  
   0.1744409294521  
   1.3572923111759  


eps = 0.700000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.7
Polinomio característico = 

                  2   3
  -1.4 +8.5x -6.7x  +x 
Raices del pol = 

   5.0814763874431  
   1.4252116344529  
   0.193311978104  
Autovalores = 

   5.0814763874431  
   0.193311978104  
   1.4252116344529  


eps = 0.800000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.8
Polinomio característico = 

                2   3
  -1.6 +9x -6.8x  +x 
Raices del pol = 

   5.0952921170896  
   1.4946092483605  
   0.2100986345499  
Autovalores = 

   5.0952921170896  
   0.2100986345499  
   1.4946092483605  


eps = 0.900000
A = 

   1.  -1.   0. 
  -2.   4.  -2. 
   0.  -1.   1.9
Polinomio característico = 

                  2   3
  -1.8 +9.5x -6.9x  +x 
Raices del pol = 

   5.1097493097029  
   1.5651862472076  
   0.2250644430895  
Autovalores = 

   5.1097493097029  
   0.2250644430895  
   1.5651862472076  


eps = 1.000000
A = 

   1.  -1.   0.
  -2.   4.  -2.
   0.  -1.   2.
Polinomio característico = 

             2   3
  -2 +10x -7x  +x 
Raices del pol = 

   5.1248854197646  
   1.6366717620673  
   0.2384428181681  
Autovalores = 

   5.1248854197646  
   0.2384428181681  
   1.6366717620673  




---------Ej 4---------
A = 

   1.  -1.
   2.  -1.
Cotas de autovalores según el teorema de Gerschorin
[0.000000, 2.000000] [-3.000000, 1.000000] 
Autovalores según scilab: 

   1.i 
  -1.i 


NO CREO PODER INSERTAR ACA UN GRAFICO EN TEXTO, SI SE TE OCURRE COMO AVISAME

---------Ej 5---------
A = 

   6.   4.   4.   1.
   4.   6.   1.   4.
   4.   1.   6.   4.
   1.   4.   4.   6.
Autovalor dominante = 15.000000
Autovector asociado = 

   1.
   1.
   0.9999999982792
   0.9999999982792
Calculado en 19 iteraciones
Autovalor dominante según scilab = 15.000000


A = 

   12.   1.   3.   4.
   1.   -3.   1.   5.
   3.    1.   6.  -2.
   4.    5.  -2.  -1.
Autovalor dominante = 14.201006
Autovector asociado = 

   1.
   0.1558606163349
   0.3183534238869
   0.2725211998613
Calculado en 31 iteraciones
Autovalor dominante según scilab = 14.201006
*/
