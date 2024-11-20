chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p7.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
// Punto a)
x = [0,0.2,0.4,0.6]
y = [1,1.2214,1.4918,1.8221]
printf("X = \n")
disp(x)
printf("y = \n")
disp(y)
printf("Polinomio interpolación lineal Lagrange\n")
p_lin_lagr = interpolLagrange(x(:,1:2), y(:,1:2))
disp(p_lin_lagr)
printf("Valor aproximado de e^(1/3) = \n")
val_p_lin_lagr = horner(p_lin_lagr, 1/3)
disp(val_p_lin_lagr)
printf("\n")
printf("Polinomio interpolación cúbica Lagrange\n")
p_cub_lagr = interpolLagrange(x, y)
disp(p_cub_lagr)
printf("Valor aproximado de e^(1/3) = \n")
val_p_cub_lagr = horner(p_cub_lagr, 1/3)
disp(val_p_cub_lagr)
printf("\n")

printf("Polinomio interpolación lineal Newton\n")
p_lin_new = interpolNewton(x(:,1:2), y(:,1:2))
disp(p_lin_new)
printf("Valor aproximado de e^(1/3) = \n")
val_p_lin_new = horner(p_lin_new, 1/3)
disp(val_p_lin_new)
printf("\n")
printf("Polinomio interpolación cúbica Newton\n")
p_cub_new = interpolNewton(x, y)
disp(p_cub_new)
printf("Valor aproximado de e^(1/3) = \n")
val_p_cub_new = horner(p_cub_new, 1/3)
disp(val_p_cub_new)

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
x = 2:0.1:2.5
y = [0.2239 0.1666 0.1104 0.0555 0.0025 -0.0484]
p = interpolNewton(x, y)
printf("Aproximación de J0(2.15) = \n")
v = horner(p, 2.15)
disp(v)
printf("Aproximación de J0(2.35) = \n")
v = horner(p, 2.35)
disp(v)

// Ejercicio 5
// En papel

// Ejercicio 6
// En papel

// Ejercicio 7
printf("\n\n---------Ej 7---------\n");
x = [0 0.15 0.31 0.5 0.6 0.75]
y = [1 1.004 1.031 1.117 1.223 1.422]
printf("x = \n")
disp(x)
printf("y = \n")
disp(y)

for n=1:3
    printf("-- Grado %d\n", n)
    A = matrizMinCuadPol(x,n)
    printf("Matriz del método")
    disp(A)
    printf("Polinomio\n")
    [p, err] = minCuadPol(A, y)
    disp(p)
    printf("Error = \n")
    disp(err)
    printf("\n")
end

// Ejercicio 8
printf("\n\n---------Ej 8---------\n");
x = [4 4.2 4.5 4.7 5.1 5.5 5.9 6.3 6.8 7.1]
y = [102.56 113.18 130.11 142.05 167.53 195.14 224.87 256.73 299.5 326.72]

for n=1:3
    printf("-- Grado %d\n", n)
    A = matrizMinCuadPol(x,n)
    printf("Matriz del método")
    disp(A)
    printf("De rango = %f\n", det(A' * A))
    printf("Polinomio\n")
    [p, err] = minCuadPol(A, y)
    pols(n) = p
    disp(p)
    printf("Error = \n")
    disp(err)
    printf("\n")
end

xx=4:0.001:7.2
plot2d(x',y,-1, leg="tabla")
plot2d(xx',[horner(pols(1),xx'),horner(pols(2),xx'),horner(pols(3),xx')],[2,3,4],leg="p1(x)@p2(x)@p3(x)")

// Ejercicio 9
printf("\n\n---------Ej 9---------\n");
printf("Ver gráfico\n")
clf()
xx = -5:0.001:5
yy = (1./(1+(xx.^2)))'
n = [2 4 6 10 14]
leyen = "f(x)"

for i=1:size(n, 2)
    grado = n(i)
    x = linspace(-5, 5, grado + 1)
    y = 1./(1+(x.^2))
    p = interpolLagrange(x, y)
    yy = [yy horner(p, xx')]
    leyen = leyen + ("@p"+string(grado)+"(x)")
end

plot2d(xx', yy, 1:size(n, 2)+1, leg=leyen)

// Ejercicio 10
printf("\n\n---------Ej 10---------\n");
clf()
x = rootsCheby(4)
// Preguntar por qué está dando una constante con Lagrange
p = interpolLagrange(x, exp(x))
p = interpolNewton(x, exp(x))
printf("Polinomio = \n")
disp(p)
xx = -1:0.001:1
plot2d(xx', (exp(xx) - horner(p, xx))')

// Ejercicio 11
printf("\n\n---------Ej 11---------\n");
printf("Ver gráfico\n")
clf()
x = rootsCheby_ab(4, 0, %pi/2)
p_cheby = interpolNewton(x, cos(x))
x = 0:%pi/(2*4):(%pi/2)
p_equidist = interpolNewton(x, cos(x))
xx = 0:0.01:(%pi/2)
plot2d(xx', [horner(p_cheby,xx)' horner(p_equidist, xx)' cos(xx)'], leg="p_cheby(x)@p_equidist(x)@cos(x)")

// Salida completa del programa
/*
---------Ej 1---------
X = 

   0.   0.2   0.4   0.6
y = 

   1.   1.2214   1.4918   1.8221
Polinomio interpolación lineal Lagrange

            
   1 +1.107x
Valor aproximado de e^(1/3) = 

   1.369

Polinomio interpolación cúbica Lagrange

                          2            3
   1 +1.0026667x +0.47625x  +0.2270833x 
Valor aproximado de e^(1/3) = 

   1.3955494

Polinomio interpolación lineal Newton

            
   1 +1.107x
Valor aproximado de e^(1/3) = 

   1.369

Polinomio interpolación cúbica Newton

                          2            3
   1 +1.0026667x +0.47625x  +0.2270833x 
Valor aproximado de e^(1/3) = 

   1.3955494


---------Ej 4---------
Aproximación de J0(2.15) = 

   0.1383688
Aproximación de J0(2.35) = 

   0.0287313


---------Ej 7---------
x = 

   0.   0.15   0.31   0.5   0.6   0.75
y = 

   1.   1.004   1.031   1.117   1.223   1.422
-- Grado 1
Matriz del método
   1.   0.  
   1.   0.15
   1.   0.31
   1.   0.5 
   1.   0.6 
   1.   0.75
Polinomio

                       
   0.929514 +0.5281021x
Error = 

   0.1567356

-- Grado 2
Matriz del método
   1.   0.     0.    
   1.   0.15   0.0225
   1.   0.31   0.0961
   1.   0.5    0.25  
   1.   0.6    0.36  
   1.   0.75   0.5625
Polinomio

                                   2
   1.011341 -0.3256988x +1.1473303x 
Error = 

   0.0307449

-- Grado 3
Matriz del método
   1.   0.     0.       0.      
   1.   0.15   0.0225   0.003375
   1.   0.31   0.0961   0.029791
   1.   0.5    0.25     0.125   
   1.   0.6    0.36     0.216   
   1.   0.75   0.5625   0.421875
Polinomio

                                   2            3
   1.0004398 -0.001541x -0.0115057x  +1.0210226x 
Error = 

   0.0105469



---------Ej 8---------
-- Grado 1
Matriz del método
   1.   4. 
   1.   4.2
   1.   4.5
   1.   4.7
   1.   5.1
   1.   5.5
   1.   5.9
   1.   6.3
   1.   6.8
   1.   7.1
De rango = 107.090000
Polinomio

                        
  -194.13824 +72.084518x
Error = 

   18.138721

-- Grado 2
Matriz del método
   1.   4.    16.  
   1.   4.2   17.64
   1.   4.5   20.25
   1.   4.7   22.09
   1.   5.1   26.01
   1.   5.5   30.25
   1.   5.9   34.81
   1.   6.3   39.69
   1.   6.8   46.24
   1.   7.1   50.41
De rango = 804.413072
Polinomio

                                    2
   1.2355604 -1.1435234x +6.6182109x 
Error = 

   0.0379857

-- Grado 3
Matriz del método
   1.   4.    16.     64.    
   1.   4.2   17.64   74.088 
   1.   4.5   20.25   91.125 
   1.   4.7   22.09   103.823
   1.   5.1   26.01   132.651
   1.   5.5   30.25   166.375
   1.   5.9   34.81   205.379
   1.   6.3   39.69   250.047
   1.   6.8   46.24   314.432
   1.   7.1   50.41   357.911
De rango = 3938.623967
Polinomio

                                    2            3
   3.4290944 -2.3792211x +6.8455778x  -0.0136746x 
Error = 

   0.0229639



---------Ej 9---------
Ver gráfico


---------Ej 10---------
Polinomio = 

                                    2            3
   0.9946153 +0.9989332x +0.5429007x  +0.1751757x 


---------Ej 11---------
Ver gráfico
*/
