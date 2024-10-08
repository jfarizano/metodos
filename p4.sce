function x = triangularSuperior(A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)

    if nA<>mA then
        error('triangularSuperior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('triangularSuperior - dimensiones incompatibles entre A y b');
        abort;
    end;

    a = [A b]; // Matriz aumentada
    n = nA;

    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = a(i, i+1:n)*x(i+1:n)
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

function x = triangularInferior(A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)

    if nA<>mA then
        error('triangularInferior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('triangularInferior - dimensiones incompatibles entre A y b');
        abort;
    end;

    a = [A b]; // Matriz aumentada
    n = nA;

    // Sustitución regresiva
    x(1) = a(1,n+1)/a(1,1);
    for i = 2:n
        sumk = a(i, 1:i-1)*x(1:i-1)
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

function [x,a,contador] = gaussElim(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  

    [nA,mA] = size(A) 
    [nb,mb] = size(b)

    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;

    a = [A b]; // Matriz aumentada

    // Eliminación progresiva
    n = nA;
    contador = 0;
    for k=1:n-1
        for i=k+1:n
            a(i, k+1:n+1) = a(i, k+1:n+1) - a(k, k+1:n+1)*a(i,k)/a(k,k);
            contador = contador + 3*(n-k+1);
            a(i, 1:k) = 0; // no hace falta para calcular la solución x
        end;
    end;

    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    contador = contador + 1;
    for i = n-1:-1:1
        sumk = sum(a(i, i+1:n)*x(i+1:n));
        x(i) = (a(i,n+1)-sumk)/a(i,i);
        contador = contador + 2*(n-i) + 2;
    end;
endfunction

function [X,A_aum] = gaussElimAmpliada(A,B)
   [nA, mA] = size(A);
   [nB, mB] = size(B);

   if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nB then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;

   // Matriz aumentada
   A_aum = [A,B];

   // Eliminación progresiva
   for i = 1:(nA-1) do
     for j = (i+1):nA do
     mji = A_aum(j,i)/A_aum(i,i) 
     A_aum(j,i)=0 
     A_aum(j,(i+1):(mA+mB)) = A_aum(j,(i+1):(mA+mB)) - mji*A_aum(i,(i+1):(mA+mB))
      end
   end

   // Sustitución regresiva
   X(nA,1:mB) = A_aum(nA,(nA+1):(nA+mB))./A_aum(nA,nA)

   for i = (nA-1):-1:1 do
      X(i,1:mB) = (A_aum(i,(mA+1):(mA+mB)) - A_aum(i,(i+1):mA)*X((i+1):mA,1:mB))./A_aum(i,i)
   end
endfunction

function detA = gaussElimDet(A)
    [nA,mA] = size(A) 

    if nA<>mA then
        error('gausselimdet - La matriz A debe ser cuadrada');
        abort;
    end;

    a = A;

    // Eliminación progresiva
    n = nA;
    for k=1:n-1
        for i=k+1:n
            a(i, k+1:n) = a(i, k+1:n) - a(k, k+1:n)*a(i,k)/a(k,k);
        end;
    end;

    // La determinante es igual a multiplicar los elementos de la diagonal
    // después de obtener la forma triangular superior
    detA = prod(diag(a))
endfunction

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
InvA = gaussElimAmpliada(A, eye(3,3))
disp(InvA)

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
A = [ 1   2  3;
      3  -2  1;
      4   2 -1]
detA = gaussElimDet(A)
printf("Determinante segun yo = %f, det segun scilab = %f\n", detA, det(A));
