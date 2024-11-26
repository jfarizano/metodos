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
        xi = a + i*h
        I = I + f(xi)
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
        xi = a + i*h
        I = I + (f(xi) * (modulo(i, 2)*2 + 2))
    end
    I = (h/3) * I
endfunction

function I = trapecioCompuestoBidim(fun, a, b, cfun, dfun, n, m)
    deff("z=f(x, y)", "z="+fun);
    deff("y=c(x)", "y="+cfun);
    deff("y=d(x)", "y="+dfun);
    h = (b - a) / n

    I = (tcb_aux(f, a, c(a), d(a), m) + tcb_aux(f, b, c(b), d(b), m)) / 2
    for i = 1:(n-1)
        xi = a + i*h
        I = I + tcb_aux(f, xi, c(xi), d(xi), m)
    end
    I = h * I
endfunction

// Igual al trapecio compuesto pero en este caso particular f
// es una función ya definida que recibe dos argumentos y el primero
// será siempre x
function I = tcb_aux(f, x, a, b, n)
    h = (b - a) / n

    I = (f(x, a) + f(x, b)) / 2
    for i = 1:(n-1)
        xi = a + i*h
        I = I + f(x, xi)
    end
    I = h * I
endfunction

function I = simpsonCompuestoBidim(fun, a, b, cfun, dfun, n, m)
    deff("z=f(x, y)", "z="+fun);
    deff("y=c(x)", "y="+cfun);
    deff("y=d(x)", "y="+dfun);
    h = (b - a) / n

    I = scb_aux(f, a, c(a), d(a), m) + scb_aux(f, b, c(b), d(b), m)
    for i = 1:(n-1)
        xi = a + i*h
        I = I + (scb_aux(f, xi, c(xi), d(xi), m) * (modulo(i, 2)*2 + 2))
    end
    I = (h/3) * I
endfunction

// Igual a simpson compuesto pero en este caso particular f
// es una función ya definida que recibe dos argumentos y el primero
// será siempre x
function I = scb_aux(f, x, a, b, n)
    h = (b - a) / n

    I = f(x, a) + f(x, b)
    for i = 1:(n-1)
        xi = a + i*h
        I = I + (f(x, xi) * (modulo(i, 2)*2 + 2))
    end
    I = (h/3) * I
endfunction
