function I = trapecio(fun, a, b)
    deff("y=f(x)", "y="+fun);
    h = b - a
 
    I = (h/2) * (f(a) + f(b))
endfunction

function I = trapecioCompuesto(fun, a, b, n)
    deff("y=f(x)", "y="+fun);
    h = (b - a) / n

    I = (f(a) + f(b)) / 2
    for i = 1:(n-1)
        x = a + i*h
        I = I + f(x)
    end
    I = h * I
endfunction

function I = simpson(fun, a, b)
    deff("y=f(x)", "y="+fun);
    h = (b - a) / 2

    I = (h/3) * (f(a) + 4*f(a+h) + f(b))
endfunction

function I = simpsonCompuesto(fun, a, b, n)
    deff("y=f(x)", "y="+fun);
    h = (b - a) / n

    I = f(a) + f(b)
    // puntos intermedios
    for i = 1:(n-1)
        x = a + i*h
        I = I + (f(x) * (modulo(i, 2)*2 + 2))
    end
    I = (h/3) * I
endfunction

