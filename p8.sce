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
