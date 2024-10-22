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
        error('gausselimampliada - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nB then
        error('gausselimampliada - dimensiones incompatibles entre A y b');
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

function [x,a] = gaussElimPP(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

    [nA,mA] = size(A) 
    [nb,mb] = size(b)

    if nA<>mA then
        error('gausselimpp - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselimpp - dimensiones incompatibles entre A y b');
        abort;
    end;

    a = [A b]; // Matriz aumentada
    n = nA;    // Tamaño de la matriz

    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i=k+1:n
            if abs(a(i,k))>amax then
                kpivot = i; amax = a(i,k);
            end;
        end;
        temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;

        for i=k+1:n
            a(i, k+1:n+1) = a(i, k+1:n+1) - a(k, k+1:n+1)*a(i,k)/a(k,k);
            a(i, 1:k) = 0; // no hace falta para calcular la solución x
        end;
    end;

    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = sum(a(i, i+1:n)*x(i+1:n));
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

function[x,a,contador]=gaussElimPPTri(A,b)
    [nA,mA]=size(A)
    [nb,mb]=size(b)
    if nA<>mA then error('gausselimpptri-La matriz A debe ser cuadrada');
                   abort;
    else 
        if mA<>nb then error('gausselimpptri - dimensiones incompatibles entre A y b');
                       abort;
        end;
    end;
    a=[A b];//Matriz aumentada
    n=nA;//Tamaño de la matriz
    contador=0;
    //Eliminación progresiva con pivoteo parcial
    for k=1:n-1 do
        kpivot=k;
        amax=abs(a(k,k)); //pivoteo
        if abs(a(k+1,k))>amax then  temp=a(k+1,:);
                                    a(k+1,:)=a(k,:);
                                    a(k,:)=temp;//Intercambiamos las filas k y k+1
        end;
        for j=k+1:n+1 do
            a(k+1,j)=a(k+1,j)-a(k,j)*a(k+1,k)/a(k,k);
            contador=contador+3;
        end;
    end;
    //Sustitución regresiva
    x(n)=a(n,n+1)/a(n,n);
    x(n-1)=(a(n-1,n+1)-a(n-1,n)*x(n))/a(n-1,n-1);
    contador=contador+4;
    for i=n-2:-1:1 do
        x(i)=(a(i,n+1)-a(i,i+1)*x(i+1)-a(i,i+2)*x(i+2))/a(i,i);
        contador=contador+5;
    end
endfunction

function [x,a,contador] = gaussElimTri(A,b)
    [nA,mA]=size(A);
    [nb,mb]=size(b);
    if nA<>mA then 
        error('gausselim-LamatrizAdebesercuadrada');
        abort;
    elseif mA<>nb then 
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
 
    a=[A b];//Matriz aumentada
    n=nA;//Tamaño de la matriz
    contador=0;
    //Eliminación progresiva sin pivoteo
    for k=1:n-1 do
        mkk=a(k+1,k)/a(k,k);
        a(k+1,k+1)=a(k+1,k+1)-mkk*a(k,k+1);
        a(k+1,n+1)=a(k+1,n+1)-mkk*a(k,n+1);
        contador=contador+5;
    end;
    //Sustitución regresiva
    x(n)=a(n,n+1)/a(n,n);
    contador=contador+1;
    for k=n-1:-1:1 do
        x(k)=(a(k,n+1)-a(k,k+1)*x(k+1))/a(k,k);
        contador=contador+3;
    end;
endfunction

function [P, L, U] = gaussElimPPLU(A)
    [nA,mA] = size(A) 

    if nA<>mA then
        error('gausselimpplu - La matriz A debe ser cuadrada');
        abort;
    end;

    a = A;
    n = nA;    // Tamaño de la matriz

    U = A
    L = eye(A)
    P = eye(A)

    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end;
        end;

        temp = U(kpivot,k:n); U(kpivot,k:n) = U(k,k:n); U(k,k:n) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:);P(k,:) = temp;

        L(k+1:n,k) = U(k+1:n,k)/U(k,k)
        U(k+1:n,k:n) = U(k+1:n,k:n) - L(k+1:n,k)*U(k, k:n) 
    end;
endfunction

function x = resolverPLU(P, L, U, b)
    y = triangularInferior(L, P * b)
    x = triangularSuperior(U, y)
endfunction

function [L, U] = doolittleLU(A)
    [nA,mA] = size(A) 

    if nA<>mA then
        error('doolittlelu - La matriz A debe ser cuadrada');
        abort;
    end;

    n = nA;    // Tamaño de la matriz

    L = eye(n,n)
    U = eye(n,n)

    U(1, :) = A(1, :)
    for i = 2:n
        L(i, 1) = A(i,1)/U(1,1)
    end

    for i=2:n
        for j = 2:n
            if i > j
                L(i, j) = (A(i, j) - L(i, 1:j-1)*U(1:j-1, j))/U(j, j)
            else
                U(i, j) = A(i, j) - (L(i, 1:i-1)*U(1:i-1, j))
            end
        end
    end
endfunction

function [U,ind] = cholesky(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('cholesky - La matriz A debe ser cuadrada');
        abort;
    end;
    // Factorización de Cholesky.
    // Trabaja únicamente con la parte triangular superior.
    //
    // ind = 1  si se obtuvo la factorización de Cholesky.
    //     = 0  si A no es definida positiva
    //
    //******************
    eps = 1.0e-8
    //******************

    n = nA
    U = zeros(n,n)

    t = A(1,1)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(1,1) = sqrt(t)
    for j = 2:n
        U(1,j) = A(1,j)/U(1,1)
    end
 
    for k = 2:n
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        if t <= eps then
            printf('Matriz no definida positiva.\n')
            ind = 0
            return
        end
        U(k,k) = sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
        end
    end
    ind = 1
endfunction

function [Q, R] = QR(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('qr - La matriz A debe ser cuadrada');
        abort;
    end;

    v(1) = norm(A(:,1))
    Q(:,1) = A(:, 1)/v(1)
    n = nA

    for k=2:n
        for j=1:k-1
            aux2(:,j) = A(:,k)'*Q(:,j)*Q(:,j)
        end
        aux = A(:,k) - sum(aux2, 2)
        v(k) = norm(aux)
        Q(:,k) =  aux / v(k)
    end

    for i=1:n
        R(i,i) = v(i)
        for j=i+1:n
            R(i,j) = A(:, j)' * Q(:, i)
        end
    end
endfunction

function [L,U] = crout(A)
    n = size(A, 1)
    L = eye(n,n)
    U = eye(n,n)
    L(:, 1) = A(:, 1)

    for j=2:n
        U(1, j) = A(1,j)/L(1,1)
    end

    for i=2:n
        for j=2:n
            if i >= j
                L(i, j) = (A(i, j) - L(i, 1:j-1)*U(1:j-1, j))
            else
                U(i, j) = (A(i, j) - (L(i, 1:i-1)*U(1:i-1, j)))/L(i, i)
            end
        end
    end
endfunction
