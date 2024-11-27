chdir(home + "/lcc/materias/metodos/lab")
exec("solo_metodos/errores.sce", -1)
exec("solo_metodos/p8.sce", -1)

// Ejercicio 1
printf("\n\n---------Ej 1---------\n");
f = "1/x"
a = 1
b = 2
printf("Integral de f(x) = %s con a = %f y b = %f\n", f, a, b)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio compuesto con n = 4:")
I_trap = trapecioCompuesto(f, a, b, 4)
disp(I_trap)
printf("Método del trapecio compuesto con n = 409:")
I_trap = trapecioCompuesto(f, a, b, 409)
disp(I_trap)
imprimir_errores(log(2), I_trap)

// Ejercicio 2
printf("\n\n---------Ej 2---------\n");
f = "%e^x * sin(x)"
a = 1
b = 3
printf("Integral de f(x) = %s con a = %f y b = %f\n", f, a, b)
printf("Según scilab:")
I_sci = integrate(f, "x", a, b)
disp(I_sci)
printf("Método del trapecio compuesto con n = 16:")
I_trap = trapecioCompuesto(f, a, b, 16)
disp(I_trap)
imprimir_errores(I_sci, I_trap)
