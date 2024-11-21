chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/p8.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
f = "log(x)"
a = 1
b = 2
printf("Integral de f(x) = %s con a = %f y b = %f\n", f, a, b)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecio(f, a, b)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpson(f, a, b)
disp(I_simp)
printf("\n")

f = "x^(1/3)"
a = 0
b = 0.1
printf("Integral de f(x) = %s con a = %f y b = %f\n", f, a, b)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecio(f, a, b)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpson(f, a, b)
disp(I_simp)
printf("\n")

f = "sin(x)^2"
a = 0
b = %pi/3
printf("Integral de f(x) = %s con a = %f y b = %f\n", f, a, b)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecio(f, a, b)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpson(f, a, b)
disp(I_simp)
printf("\n")

// Ejercicios 2 y 3
printf("\n\n---------Ejs 2 y 3---------\n");
f = "1/x"
a = 1
b = 3
n = 4
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")

f = "x^3"
a = 0
b = 2
n = 4
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")

f = "x*((1+x^2)^(1/2))"
a = 0
b = 3
n = 6
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")

f = "sin(%pi*x)"
a = 0
b = 1
n = 8
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")

f = "x * sin(x)"
a = 0
b = 2*%pi
n = 8
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)

printf("\n")
f = "x^2 * exp(x)"
a = 0
b = 1
n = 8
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")

// Ejercicio 4
printf("\n\n---------Ej 4---------\n");
f = "1/(x+1)"
a = 0
b = 1.5
n = 10
I_real = 0.9162907
printf("Integral de f(x) = %s con a = %f, b = %f, n = %d\n", f, a, b, n)
printf("Valor real:")
disp(I_real)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio:")
I_trap = trapecioCompuesto(f, a, b, n)
disp(I_trap)
printf("Método de Simpson:")
I_simp = simpsonCompuesto(f, a, b, n)
disp(I_simp)
printf("\n")


// Salida completa del programa
/*
---------Ej 1---------
Integral de f(x) = log(x) con a = 1.000000 y b = 2.000000
Según scilab:
   0.3862944
Método del trapecio:
   0.3465736
Método de Simpson:
   0.3858346

Integral de f(x) = x^(1/3) con a = 0.000000 y b = 0.100000
Según scilab:
   0.0348119
Método del trapecio:
   0.0232079
Método de Simpson:
   0.0322962

Integral de f(x) = sin(x)^2 con a = 0.000000 y b = 1.047198
Según scilab:
   0.3070924
Método del trapecio:
   0.3926991
Método de Simpson:
   0.3054326



---------Ejs 2 y 3---------
Integral de f(x) = 1/x con a = 1.000000, b = 3.000000, n = 4
Según scilab:
   1.0986123
Método del trapecio:
   1.1166667
Método de Simpson:
   1.1

Integral de f(x) = x^3 con a = 0.000000, b = 2.000000, n = 4
Según scilab:
   4.
Método del trapecio:
   4.25
Método de Simpson:
   4.

Integral de f(x) = x*((1+x^2)^(1/2)) con a = 0.000000, b = 3.000000, n = 6
Según scilab:
   10.207592
Método del trapecio:
   10.312201
Método de Simpson:
   10.206346

Integral de f(x) = sin(%pi*x) con a = 0.000000, b = 1.000000, n = 8
Según scilab:
   0.6366198
Método del trapecio:
   0.6284174
Método de Simpson:
   0.6367055

Integral de f(x) = x * sin(x) con a = 0.000000, b = 6.283185, n = 8
Según scilab:
  -6.2831853
Método del trapecio:
  -5.9568332
Método de Simpson:
  -6.2975102

Integral de f(x) = x^2 * exp(x) con a = 0.000000, b = 1.000000, n = 8
Según scilab:
   0.7182818
Método del trapecio:
   0.7288902
Método de Simpson:
   0.7183215



---------Ej 4---------
Integral de f(x) = 1/(x+1) con a = 0.000000, b = 1.500000, n = 10
Valor real:
   0.9162907
Según scilab:
   0.9162907
Método del trapecio:
   0.9178617
Método de Simpson:
   0.9163064
*/
