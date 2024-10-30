chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p4.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
A = [ 3 -2 -1; 
      0 -2  2; 
      0  0  1]
b = [0 6 -1]'
x = triangularSuperior(A,b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)

printf("\n\n")

A = [ 1  0 0; 
      2 -2 0; 
     -1 -2 3]
b = [-1 6 0]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
x = triangularInferior(A, b)
disp(x)

// Ejercicio 2
printf("\n\n---------Ej 2---------\n");
A = [ 1  1  0  3;
      2  1 -1  1;
      3 -1 -1  2;
     -1  2  3 -1]
b = [4 1 -3 4]'
x_real=[-1 2 0 1]'
x = gaussElim(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("x real = \n")
disp(x_real)
printf("\n\n")

A = [ 1 -1  2 -1;
      2 -2  3 -3;
      1  1  1  0;
      1 -1  4  3]
b = [-8 -20 -2 -4]'
x = gaussElim(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")

A = [ 1  1  0  4;
      2  1 -1  1;
      4 -1 -2  2;
      3 -1 -1  2]
b = [2 1 0 -3]'
x = gaussElim(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")

// Ejercicio 3
printf("\n\n---------Ej 3---------\n");
A = [ 1   2  3;
      3  -2  1;
      4   2 -1]
b = [ 14  9 -2;
       2 -5  2;
       5 19 12]
x = gaussElimAmpliada(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")
printf("x1 = \n")
disp(gaussElim(A, b(:,1)))
printf("x2 = \n")
disp(gaussElim(A, b(:,2)))
printf("x3 = \n")
disp(gaussElim(A, b(:,3)))
printf("Inversa de A = \n")
InvA = gaussElimAmpliada(A, eye(A))
disp(InvA)

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
A = [ 1   2  3;
      3  -2  1;
      4   2 -1]
detA = gaussElimDet(A)
printf("Determinante segun yo = %f, det segun scilab = %f\n", detA, det(A));

// Ejercicio 5
printf("\n\n---------Ej 5---------\n");
A = [0 2 3; 2 0 3; 8 16 -1]
b = [7 13 -3]'
[x,a] = gaussElimPP(A,b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")

A = [ 1  1  0  3;
      2  1 -1  1;
      3 -1 -1  2;
     -1  2  3 -1]
b = [4 1 -3 4]'
x_real=[-1 2 0 1]'
x = gaussElimPP(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("x real = \n")
disp(x_real)
printf("\n\n")

A = [ 1 -1  2 -1;
      2 -2  3 -3;
      1  1  1  0;
      1 -1  4  3]
b = [-8 -20 -2 -4]'
x = gaussElimPP(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")

A = [ 1  1  0  4;
      2  1 -1  1;
      4 -1 -2  2;
      3 -1 -1  2]
b = [2 1 0 -3]'
x = gaussElimPP(A, b)
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
printf("x = \n")
disp(x)
printf("\n\n")

// Ejercicio 6
printf("\n\n---------Ej 6---------\n");
A = [1 1 0 0;
     1 1 1 0;
     0 1 1 1; 
     0 0 1 1]
b = [2 3 3 2]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
x = gaussElimTri(A, b)
printf("x sin pp = \n")
disp(x)
x = gaussElimPPTri(A, b)
printf("x con pp = \n")
disp(x)
printf("\n\n")

A = [2 2 0 0; 
     5 2 2 0; 
     0 5 2 2; 
     0 0 5 2]
b = [6 15 24 23]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
x = gaussElimTri(A, b)
printf("x sin pp = \n")
disp(x)
x = gaussElimPPTri(A, b)
printf("x con pp = \n")
disp(x)
printf("\n\n")

A = [2 1 0 0;
     1 2 1 0;
     0 1 2 1;
     0 0 1 2]
b = [1 1 1 1]'
printf("A = \n")
disp(A)
printf("b = \n")
disp(b)
x = gaussElimTri(A, b)
printf("x sin pp = \n")
disp(x)
x = gaussElimPPTri(A, b)
printf("x con pp = \n")
disp(x)
printf("\n\n")

// Ejercicio 7
printf("\n\n---------Ej 7---------\n");
A = [2 1 1 0;
     4 3 3 1;
     8 7 9 5;
     6 7 9 8]
printf("A = \n")
disp(A)
[P, L, U] = gaussElimPPLU(A)
printf("P = \n")
disp(P)
printf("L = \n")
disp(L)
printf("U = \n")
disp(U)

// Ejercicio 8
printf("\n\n---------Ej 8---------\n");
A = [ 1.012 -2.132  3.104;
     -2.132  4.096 -7.013;
      3.104 -7.013  0.014]
printf("A = \n")
disp(A)
[P, L, U] = gaussElimPPLU(A)
printf("P propia = \n")
disp(P)
printf("L propia = \n")
disp(L)
printf("U propia = \n")
disp(U)
[L, U, P] = lu(A)
printf("P scilab = \n")
disp(P)
printf("L scilab = \n")
disp(L)
printf("U scilab = \n")
disp(U)

A = [ 2.1756 4.0231 -2.1732 5.1967
     -4.0231 6.0000  0.0000 1.1973;
     -1.0000 5.2107  1.1111 0.0000;
      6.0235 7.0000  0.0000 4.1561]
printf("A = \n")
disp(A)
[P, L, U] = gaussElimPPLU(A)
printf("P propia = \n")
disp(P)
printf("L propia = \n")
disp(L)
printf("U propia = \n")
disp(U)
[L, U, P] = lu(A)
printf("P scilab = \n")
disp(P)
printf("L scilab = \n")
disp(L)
printf("U scilab = \n")
disp(U)

// Ejercicio 9
printf("\n\n---------Ej 9---------\n");
A = [1  2  -2   1;
     4  5  -7   6;
     5  25 -15 -3;
     6 -12 -6   22]
b = [1 2 0 1]'
printf("A = \n")
disp(A)
[P, L, U] = gaussElimPPLU(A)
printf("P = \n")
disp(P)
printf("L = \n")
disp(L)
printf("U = \n")
disp(U)
printf("b = \n")
disp(b)
x = resolverPLU(P, L, U, b)
printf("x = \n")
disp(x)
b = [2 2 1 0]'
printf("nueva b = \n")
disp(b)
x = resolverPLU(P, L, U, b)
printf("x = \n")
disp(x)

// Ejercicio 10
printf("\n\n---------Ej 10---------\n");
A = [1 2  3  4;
     1 4  9  16;
     1 8  27 64;
     1 16 81 256]
b = [2 10 44 190]'
printf("A = \n")
disp(A)
[L, U] = doolittleLU(A)
printf("L = \n")
disp(L)
printf("U = \n")
disp(U)
printf("b = \n")
disp(b)
x = resolverPLU(eye(A), L, U, b)
printf("x = \n")
disp(x)

// Ejercicio 11
printf("\n\n---------Ej 11---------\n");
A = [16 -12  8 -16;
     12  18 -6  9;
     8  -6   5 -10;
     16  9  -10 46]
printf("A = \n")
disp(A)
[UA, ind] = cholesky(A)
if (ind == 0) then
    printf("No se pudo obtener la factorización de Cholesky\n")
else
    printf("U = \n")
    disp(UA)
    if (UA == chol(A)) then
        printf("La factorización es correcta\n")
    else
        printf("La factorización es incorrecta\n")
    end
end

B = [4 1 1;
     8 2 2;
     1 2 3]
printf("B = \n")
disp(B)
[UB, ind] = cholesky(B)
if (ind == 0) then
    printf("No se pudo obtener la factorización de Cholesky\n")
else
    printf("U = \n")
    disp(UB)
    if (UB == chol(B)) then
        printf("La factorización es correcta\n")
    else
        printf("La factorización es incorrecta\n")
    end
end

C = [1 2;
     2 4]
printf("C = \n")
disp(C)
[U, ind] = cholesky(C)
if (ind == 0) then
    printf("No se pudo obtener la factorización de Cholesky\n")
else
    printf("C = \n")
    disp(C)
    if (UC == chol(C)) then
        printf("La factorización es correcta\n")
    else
        printf("La factorización es incorrecta\n")
    end
end


// Ejercicio 12
// item a)
printf("\n\n---------Ej 12---------\n");
A = [16 -12  8;
    -16  18 -6;
     8  -6   8]
printf("A = \n")
disp(A)
b = [76 -66 46]'
[U, ind] = cholesky(A)
if (ind == 0) then
    printf("No se pudo obtener la factorización de Cholesky\n")
else
    printf("U = \n")
    disp(U)
    printf("b = \n")
    disp(b)
    x = resolverPLU(eye(A), U', U)
    printf("x = \n")
    disp(x)
end

// item b)
[Q, R] = QR(A)
printf("Q = \n")
disp(Q)
printf("R = \n")
disp(R)
printf("b = \n")
disp(b)
x = triangularSuperior(R, Q' * b)
printf("x = \n")
disp(x)


// Salida completa del programa
/*
---------Ej 1---------
A = 

   3.  -2.  -1.
   0.  -2.   2.
   0.   0.   1.
b = 

   0.
   6.
  -1.
x = 

  -3.
  -4.
  -1.


A = 

   1.   0.   0.
   2.  -2.   0.
  -1.  -2.   3.
b = 

  -1.
   6.
   0.
x = 

  -1.
  -4.
  -3.


---------Ej 2---------
A = 

   1.   1.   0.   3.
   2.   1.  -1.   1.
   3.  -1.  -1.   2.
  -1.   2.   3.  -1.
b = 

   4.
   1.
  -3.
   4.
x = 

  -1.
   2.
   0.
   1.
x real = 

  -1.
   2.
   0.
   1.


A = 

   1.  -1.   2.  -1.
   2.  -2.   3.  -3.
   1.   1.   1.   0.
   1.  -1.   4.   3.
b = 

  -8.
  -20.
  -2.
  -4.
x = 

   Nan
   Nan
   Nan
   Nan


A = 

   1.   1.   0.   4.
   2.   1.  -1.   1.
   4.  -1.  -2.   2.
   3.  -1.  -1.   2.
b = 

   2.
   1.
   0.
  -3.
x = 

  -4.
   0.6666666666667
  -7.
   1.3333333333333




---------Ej 3---------
A = 

   1.   2.   3.
   3.  -2.   1.
   4.   2.  -1.
b = 

   14.   9.   -2. 
   2.   -5.    2. 
   5.    19.   12.
x = 

   1.   2.   2.
   2.   5.   1.
   3.  -1.  -2.


x1 = 

   1.
   2.
   3.
x2 = 

   2.
   5.
  -1.
x3 = 

   2.
   1.
  -2.
Inversa de A = 

   0.      0.1428571428571   0.1428571428571
   0.125  -0.2321428571429   0.1428571428571
   0.25    0.1071428571429  -0.1428571428571


---------Ej 4---------
Determinante segun yo = 56.000000, det segun scilab = 56.000000


---------Ej 5---------
A = 

   0.   2.    3.
   2.   0.    3.
   8.   16.  -1.
b = 

   7.
   13.
  -3.
x = 

   2.
  -1.
   3.


A = 

   1.   1.   0.   3.
   2.   1.  -1.   1.
   3.  -1.  -1.   2.
  -1.   2.   3.  -1.
b = 

   4.
   1.
  -3.
   4.
x = 

  -1.
   2.
   1.850371708D-17
   1.
x real = 

  -1.
   2.
   0.
   1.


A = 

   1.  -1.   2.  -1.
   2.  -2.   3.  -3.
   1.   1.   1.   0.
   1.  -1.   4.   3.
b = 

  -8.
  -20.
  -2.
  -4.
x = 

  -15.
   7.
   6.
  -2.


A = 

   1.   1.   0.   4.
   2.   1.  -1.   1.
   4.  -1.  -2.   2.
   3.  -1.  -1.   2.
b = 

   2.
   1.
   0.
  -3.
x = 

  -4.
   0.6666666666667
  -7.
   1.3333333333333




---------Ej 6---------
A = 

   1.   1.   0.   0.
   1.   1.   1.   0.
   0.   1.   1.   1.
   0.   0.   1.   1.
b = 

   2.
   3.
   3.
   2.
x sin pp = 

   Nan
   Nan
   Nan
   2.
x con pp = 

   1.
   1.
   1.
   1.


A = 

   2.   2.   0.   0.
   5.   2.   2.   0.
   0.   5.   2.   2.
   0.   0.   5.   2.
b = 

   6.
   15.
   24.
   23.
x sin pp = 

   1.
   2.
   3.
   4.
x con pp = 

   1.
   2.
   3.
   4.


A = 

   2.   1.   0.   0.
   1.   2.   1.   0.
   0.   1.   2.   1.
   0.   0.   1.   2.
b = 

   1.
   1.
   1.
   1.
x sin pp = 

   0.4
   0.2
   0.2
   0.4
x con pp = 

   0.4
   0.2
   0.2
   0.4




---------Ej 7---------
A = 

   2.   1.   1.   0.
   4.   3.   3.   1.
   8.   7.   9.   5.
   6.   7.   9.   8.
P = 

   0.   0.   1.   0.
   0.   0.   0.   1.
   0.   1.   0.   0.
   1.   0.   0.   0.
L = 

   1.     0.                0.                0.
   0.75   1.                0.                0.
   0.5   -0.2857142857143   1.                0.
   0.25  -0.4285714285714   0.3333333333333   1.
U = 

   8.   7.     9.                5.             
   0.   1.75   2.25              4.25           
   0.   0.    -0.8571428571429  -0.2857142857143
   0.   0.     0.                0.6666666666667


---------Ej 8---------
A = 

   1.012  -2.132   3.104
  -2.132   4.096  -7.013
   3.104  -7.013   0.014
P propia = 

   0.   0.   1.
   0.   1.   0.
   1.   0.   0.
L propia = 

   1.                0.                0.
  -0.6868556701031   1.                0.
   0.3260309278351  -0.2142472825164   1.
U propia = 

   3.104  -7.013            0.014          
   0.     -0.720918814433  -7.0033840206186
   0.      0.               1.598979572174 
P scilab = 

   0.   0.   1.
   0.   1.   0.
   1.   0.   0.
L scilab = 

   1.                0.                0.
  -0.6868556701031   1.                0.
   0.3260309278351  -0.2142472825164   1.
U scilab = 

   3.104  -7.013            0.014          
   0.     -0.720918814433  -7.0033840206186
   0.      0.               1.598979572174 
A = 

   2.1756   4.0231  -2.1732   5.1967
  -4.0231   6.       0.       1.1973
  -1.       5.2107   1.1111   0.    
   6.0235   7.       0.       4.1561
P propia = 

   0.   0.   0.   1.
   0.   1.   0.   0.
   1.   0.   0.   0.
   0.   0.   1.   0.
L propia = 


         column 1 to 3

   1.                0.                0.             
  -0.6679007221715   1.                0.             
   0.3611853573504   0.1400243356811   1.             
  -0.1660164356271   0.596967957022   -0.5112736977729

         column 4

   0.
   0.
   0.
   1.
U propia = 

   6.0235   7.              0.       4.1561         
   0.       10.6753050552   0.       3.973162191417 
   0.       0.             -2.1732   3.1392381399097
   0.       0.              0.      -0.0768597162361
P scilab = 

   0.   0.   0.   1.
   0.   1.   0.   0.
   1.   0.   0.   0.
   0.   0.   1.   0.
L scilab = 


         column 1 to 3

   1.                0.                0.             
  -0.6679007221715   1.                0.             
   0.3611853573504   0.1400243356811   1.             
  -0.1660164356271   0.596967957022   -0.5112736977729

         column 4

   0.
   0.
   0.
   1.
U scilab = 

   6.0235   7.              0.       4.1561         
   0.       10.6753050552   0.       3.973162191417 
   0.       0.             -2.1732   3.1392381399097
   0.       0.              0.      -0.0768597162361


---------Ej 9---------
A = 

   1.   2.   -2.    1. 
   4.   5.   -7.    6. 
   5.   25.  -15.  -3. 
   6.  -12.  -6.    22.
P = 

   0.   0.   0.   1.
   0.   0.   1.   0.
   0.   1.   0.   0.
   1.   0.   0.   0.
L = 

   1.                0.                0.    0.
   0.8333333333333   1.                0.    0.
   0.6666666666667   0.3714285714286   1.    0.
   0.1666666666667   0.1142857142857   0.2   1.
U = 

   6.  -12.  -6.                22.            
   0.   35.  -10.              -21.333333333333
   0.   0.    0.7142857142857  -0.7428571428571
   0.   0.    0.               -0.08           
b = 

   1.
   2.
   0.
   1.
x = 

   9.8333333333334
  -6.1666666666667
  -5.5
  -7.5
nueva b = 

   2.
   2.
   1.
   0.
x = 

   19.5
  -17.
  -18.
  -19.5


---------Ej 10---------
A = 

   1.   2.    3.    4.  
   1.   4.    9.    16. 
   1.   8.    27.   64. 
   1.   16.   81.   256.
P = 

   0.   0.   0.   1.
   0.   0.   1.   0.
   0.   1.   0.   0.
   1.   0.   0.   0.
L = 

   1.   0.   0.   0.
   1.   1.   0.   0.
   1.   3.   1.   0.
   1.   7.   6.   1.
U = 

   1.   2.   3.   4. 
   0.   2.   6.   12.
   0.   0.   6.   24.
   0.   0.   0.   24.
b = 

   2.
   10.
   44.
   190.
x = 

  -1.
   1.
  -1.
   1.


---------Ej 11---------
A = 

   16.  -12.   8.   -16.
   12.   18.  -6.    9. 
   8.   -6.    5.   -10.
   16.   9.   -10.   46.
U = 

   4.  -3.   2.  -4.
   0.   3.   0.  -1.
   0.   0.   1.  -2.
   0.   0.   0.   5.
La factorización es correcta
B = 

   4.   1.   1.
   8.   2.   2.
   1.   2.   3.
U = 

   2.   0.5               0.5            
   0.   1.3228756555323   1.3228756555323
   0.   0.                1.             
La factorización es correcta
C = 

   1.   2.
   2.   4.
Matriz no definida positiva.
No se pudo obtener la factorización de Cholesky


---------Ej 12---------
A = 

   16.  -12.   8.
  -16.   18.  -6.
   8.   -6.    8.
U = 

   4.  -3.   2.
   0.   3.   0.
   0.   0.   2.
b = 

   76.
  -66.
   46.
x = 

   3.
  -1.
   2.
Q = 

   0.6666666666667   0.5962847939999  -0.4472135955   
  -0.6666666666667   0.7453559924999   0.             
   0.3333333333333   0.298142397       0.8944271909999
R = 

   24.  -22.               12.            
   0.    4.4721359549996   2.6832815729997
   0.    0.                3.5777087639997
b = 

   76.
  -66.
   46.
x = 

   4.5
   1.
   2.
*/
