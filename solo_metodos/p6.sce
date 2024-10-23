function [valor, z1, contador] = potencia(A, z0, tol, iter)
    contador = 0
    delta = %inf

    while delta > tol && contador < iter
        w = A*z0
        z1 = w / norm(w, %inf)
        [m, j] = max(abs(w))
        valor = w(j) / z0(j)
        z1 = w / valor

        delta = norm(z0 - z1, %inf)
        z0 = z1
        contador = contador + 1
    end
endfunction

function valor = maxAutovalor(A)
    valores = eigs(A)
    [m, i] = max(abs(valores))
    valor = valores(i)
end

function circ = circulosGers(A)
    // Dada una matrix nxn, devuelve una matrix n x 2
    // Donde cada fila es (c, r) que corresponden al centro y radio de un círculo
    circ(:,1) = diag(A)
    circ(:,2) = sum(abs(A - diag(diag(A))), 2)'
endfunction

function dibujarCirculosGers(A)
    circ = circulosGers(A)
    n = size(A, 1)
    max_y = ceil(max(circ(:, 2)))
    aux1 = circ(:, 1) + circ(:, 2)
    aux2 = circ(:, 1) - circ(:, 2)
    max_x = ceil(max(max(abs(aux1)), max(abs(aux2))))
    bordes = [-max_x -max_y/2 max_x max_y/2]
    plot2d(0, 0, 0, "031", " ", bordes, frameflag=3, axesflag=5)
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";

    for i=1:n
        // (x, y) representa la esquina sup. izq. del rectangulo que contiene al circ
        r = circ(i, 2)
        c = circ(i, 1)
        x = c - r
        y = r
        xarc(x, y, r*2, r*2, 0, 360*64)
    end
endfunction

function dibujarCirculosGersValor(A)
    dibujarCirculosGers(A)
    valores = spec(A)
    xpoly(real(valores), imag(valores), "marks")
    e=gce()
    set(e,"mark_style",3)
endfunction
