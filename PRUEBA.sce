// función que te calcula el Jacobiano según la posición del end-efector
function [J] = jacobianos(matrizJ2,matrizJ3,i)
    s1 = [0; 0; 1];
    s2 = [matrizJ2(i,2); -matrizJ2(i,1); 1];
    s3 = [matrizJ3(i,2); -matrizJ3(i,1); 1];
    J = [s1,s2,s3]
endfunction

// función fordward. No la utilizamos para nuestro ejercicio
/*function [T] = forward(J,g)
    T = J*g
endfunction*/

// función que calcula el determinante y te dice si es invertible o no para evitar SINGULARIDADES
function [detJ] = determinante(J)
    detJ = det(J);
    if detJ == 0 then
        printf('No invertible\n');
    else
        printf('Invertible\n');
    end
endfunction

// función para invertir el Jacobiano
function [invJ] = invertida(J)
    invJ = inv(J)
endfunction

// función para calcular el metodo IIKP
function [G] = inverse(invJ,T)
    G = invJ*T
endfunction


//FUNCIÓN PARA CALCULAR LA POSICIÓN DEL END-EFFECTOR SEGÚN UNA VELOCIDAD VARIABLE Y ALMACENAR DATOS EN UNA MATRIZ

function [matrizP,n] = posicion(v)
    i = 1;
    Px = 0.9;
    Py = -0.2;
    if v == 0 then
        n = 1;
    else
        n = floor(((0.7-0.2)/v)+1)+1; //nos sirve para marcar la longitud en filas que tendrá la matriz. La función floor nos dará la parte entera de la operación
    end
    matrizP = zeros(n,2); //creamos una matriz nx2 llena de 0 para luego rellenarla a nuestro gusto
    t = 0;
    if v == 0 then
        matrizP = [Px,Py];
    else
        while Py > (-0.7)
            Py = -0.2 - v*t;
            matrizP(i,1) = Px; //rellenamos la primera columna de nuestra matriz con los valores de Px que siempre es constante (0.9)
            matrizP(i,2) = Py; //rellenamos la segunda columna de nuestra matriz con los valores de Py que varia respecto el tiempo
            if Py < -0.7 then
                Py = -0.7; //limitamos el valor de Py a -0.7m que es nuestro máximo por si nos pasaramos por la velocidad del movimiento
                matrizP(i,2) = Py;
            end
            i = i+1;
            t = t+1;
        end
    end
endfunction


//FUNCIÓN QUE CALCULA LA POSICION DE LA JOINT 3 Y ALMACENA TODAS LAS POSICIONES EN UNA MATRIZ

function [matrizJ3] = posjoint3(matrizP)
    matrizJ3 = zeros(n,2);
    i = 1;
    J3x = 0.5;
    while i <= n
        matrizJ3(i,1) = J3x;
        matrizJ3(i,2) = matrizP(i,2)+0.2;
        i = i+1;
    end
endfunction


//FUNCIÓN QUE CALCULA LA POSICION DE LA JOINT 2 Y ALMACENA TODAS LAS POSICIONES EN UNA MATRIZ

function [matrizJ2] = posjoint2(matrizJ3)
    matrizJ2 = zeros(n,2);
    i = 1;
    while i <= n
        a = sqrt(matrizJ3(i,1)^2+matrizJ3(i,2)^2); //usando pitagoras calcularemos la hipotenusa del triangulo formado entre la joint 2 y 3
        x = 0.62^2-0.57^2+a^2;
        y = sqrt(0.62^2-x^2);
        teta1 = atan(y/x);
        teta2 = atan(matrizJ3(i,1)/matrizJ3(i,2));
        teta3 = teta2-teta1;
        matrizJ2(i,1) = 0.62*sin(teta3);
        matrizJ2(i,2) = 0.62*cos(teta3);
        i = i+1;
    end
endfunction


//FUNCIÓN PARA CALCULAR TODAS LA VELOCIDADES DE CADA JOINT POR CADA POSICIÓN DEL END-EFECTOR
// usamos todas las funciones anteriores para llamarlas en una sola y hacer el calculo total y graficarlo
function [P,J3,J2,W] = calculo(v)
    [matrizP,n] = posicion(v);
    P = matrizP;
    [matrizJ3] = posjoint3(matrizP);
    J3 = matrizJ3;
    [matrizJ2] = posjoint2(matrizJ3);
    J2 = matrizJ2;
    i = 1;
    W = zeros(3,n);
    while i <= n
        [J] = jacobianos(matrizJ2,matrizJ3,i);
        [detJ] = determinante(J);
        [invJ] = invertida(J);
        T = [matrizP(i,2);-matrizP(i,1);0];
        [G] = inverse(invJ,T);
        W(1,i) =  G(1);
        W(2,i) =  G(2);
        W(3,i) =  G(3);
        i = i+1;
    end
    
    // en esta parte graficaremos las velocidades de cada joint de forma variable con el tiempo
    x = W(1,:);
    y = W(2,:);
    z = W(3,:);
    time = [0:(n-1)];
    xlabel('tiempo (s)');
    ylabel('velocidades (rad/s)');
    title('Grafica de velocidades de las joints');
    for i=1:n
        set(gca(),"auto_clear","off");
        xpause(500000); // nos esperamos 500000 microsegundos en cada punto
        plot(time(1:i),x(1:i),time(1:i),y(1:i),time(1:i),z(1:i),'linewidth',2);
    end
    legend('W1','W2','W3');
endfunction


