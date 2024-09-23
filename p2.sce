function r = raicesRobustas(pol)
    a = coeff(pol, 2);
    b = coeff(pol, 1);
    c = coeff(pol, 0);

    disc = b**2 - 4 * a * c;

    if (disc <= 0) || (a == 0) then
        return %nan;
    end

    sqrt_disc = sqrt(disc);

    if (b < 0) then
        r(1) = (-b + sqrt_disc) / 2 * a;
        r(2) = 2 * c / (-b + sqrt_disc);
    else
        r(1) = 2 * c / (-b - sqrt_disc);
        r(2) = (-b - sqrt_disc) / 2 * a;
    end
endfunction

// funcion f es la ley de la función dada por un string, usa como 
// variable x
// v es el valor donde se evaluará la derivada
// n es el orden de derivación
// h es el paso de derivación

function valor = derivada(f,v,n,h)
    deff("y=DF0(x)","y="+f);
    if n==0 then valor = DF0(v);
    else
        for i=1:(n-1)
        deff("y=DF"+string(i)+"(x)","y=(DF"+string(i-1)+"(x+"+string(h)+")-DF"+string(i-1)+"(x))/"+string(h));
        end
        deff("y=DFn(x)","y=(DF"+string(n-1)+"(x+"+string(h)+")-DF"+string(n-1)+"(x))/"+string(h));
        valor = DFn(v);
    end
endfunction

// usando numderivative
// esta función utiliza un orden para numderivative igual a 4
function valor = derivadaNum(f,v,n,h)
    deff("y=DF0(x)","y="+f);
    if n==0 then valor = DF0(v);
    else
        for i=1:(n-1)
        deff("y=DF"+string(i)+"(x)","y=numderivative(DF"+string(i-1)+",x,"+string(h)+",4)");
        end
        deff("y=DFn(x)","y=numderivative(DF"+string(n-1)+",x,"+string(h)+",4)");
        valor = DFn(v);
    end
endfunction

function v = deriv(f,n,vx,h)
    if n==0 then v=f(vx)
    else v=(deriv(f,n-1,vx+h,h)-deriv(f,n-1,vx,h))/h
    end
endfunction


function [b, d] = calcularhorner(pol, x)
    n = degree(pol);
    b = coeff(pol, n);
    d = b;

    for i = n-1:-1:0
        a = coeff(pol, i);
        b =  a + x * b;

        if i > 0 then
            d = d * x + b;
        end
    end
endfunction


function y = taylor(f, a, n, h, x)
    fact = 1;
    xma = x - a;
    resta = 1;
    deff("y=fun(x)","y="+f);
    y = fun(a);
    for i = 1:n
        resta = resta * xma;
        fact = fact * i;
        y = y + resta * (1/fact) * derivadaNum(f, a, i, h);
    end
endfunction

// Ej 1
printf("\n\n---------Ej 1---------\n");
pol = poly([-0.0001 10000.0 0.0001],"x","coeff");
r = raicesRobustas(pol);
printf("Polinomio:\n");
disp(pol);
printf("Raíz positiva: \n");
disp(r(1));
