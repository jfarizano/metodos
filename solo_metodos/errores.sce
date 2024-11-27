function e = error_abs(x_real, x_aprox)
    e = abs(x_real - x_aprox)
endfunction

function e = error_rel(x_real, x_aprox)
    e = abs(x_real - x_aprox) / abs(x_real)
endfunction

function imprimir_errores(x_real, x_aprox)
    printf("Valor real = %f ~ Valor aproximado = %f\n", x_real, x_aprox)
    printf("Error absoluto =")
    disp(error_abs(x_real, x_aprox))
    printf("Error relativo =")
    disp(error_rel(x_real, x_aprox))
endfunction
