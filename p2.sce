chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p2.sce", -1)

// Ej 1
printf("\n\n---------Ej 1---------\n");
pol = poly([-0.0001 10000.0 0.0001],"x","coeff");
r = raicesRobustas(pol);
printf("Polinomio:\n");
disp(pol);
printf("Raíz positiva: \n");
disp(r(1));
printf("Raíz negativa: \n");
disp(r(2));
printf("Raices con roots:\n")
disp(roots(pol));

printf("\nOtro ejemplo: \n");
pol = poly([2 6 3],"x","coeff");
r = raicesRobustas(pol);
printf("Polinomio:\n");
disp(pol);
printf("Raíz positiva: \n");
disp(r(1));
printf("Raíz negativa: \n");
disp(r(2));
printf("Raices con roots:\n")
disp(roots(pol));

// Ej 3
printf("\n\n---------Ej 3---------\n");
pol = poly([2 6 3 7 4],"x","coeff");
printf("Polinomio de prueba:\n");
disp(pol);
x = 42;
[b, d] = calcularHorner(pol, x);
printf("Resultado de evaluar en x = %f,  p(x) = %f\n", x, b);
printf("Derivada segun horner = %f\n", d);
printf("Polinomio derivada de p:\n");
pol = poly([6 6 21 16],"x","coeff");
disp(pol);
[b, d] = calcularHorner(pol, x);
printf("Resultado de evaluar polinomio derivada = %f", b);

// Ej 6
printf("\n\n---------Ej 6---------\n");
// Funcion que calcula e^x mediante un polinomio de taylor alrededor de 0
function y = taylor_e_cero(x, n)
    fact = 1;
    xn = 1;
    y = 1;
    for i = 1:n
        xn = xn * x;
        fact = fact * i;
        y = y + xn * (1/fact);
    end
endfunction

printf("Valor real de e^(-2):\n");
disp(%e^(-2));
y = taylor_e_cero(-2, 10);
printf("Calculando a partir de e^(-2):\n");
disp(y);
y = 1 / taylor_e_cero(2, 10);
printf("Calculando a partir de 1/e^2:\n");
disp(y);


// Salida completa del programa
/*
---------Ej 1---------
Polinomio:

                          2
  -0.0001 +10000x +0.0001x 
Raíz positiva: 

   0.00000001
Raíz negativa: 

  -100000000.
Raices con roots:

  -100000000.
   0.00000001

Otro ejemplo: 
Polinomio:

            2
   2 +6x +3x 
Raíz positiva: 

  -0.4226497308104
Raíz negativa: 

  -1.5773502691896
Raices con roots:

  -1.5773502691896  
  -0.4226497308104  


---------Ej 3---------
Polinomio de prueba:

            2    3    4
   2 +6x +3x  +7x  +4x 
Resultado de evaluar en x = 42.000000,  p(x) = 12970946.000000
Derivada segun horner = 1222710.000000
Polinomio derivada de p:

             2     3
   6 +6x +21x  +16x 
Resultado de evaluar polinomio derivada = 1222710.000000

---------Ej 6---------
Valor real de e^(-2):

   0.1353352832366
Calculando a partir de e^(-2):

   0.1353791887125
Calculando a partir de 1/e^2:

   0.1353364076419

*/
