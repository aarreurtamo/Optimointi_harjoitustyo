
% Funktio, joka ratkaisee epälineaarisen optiomointiongelman newton
% gradientti menetelmää ja sakkomenetelmää hyödyntäen.
% 
% Syötteet
% 1. Optimoitava funktio f numeerisena.
% 2. Epäyhtälörajoitteet vektorimuodossa g ja muodossa g <= 0
% 3. Yhtälörajoitteet vektorimuodossa h ja muodossa h = 0
% 4. xk alkuarvaus pystyvektorina.
% 5. toleranssi tol, jota pienemmäksi gradientin arvo minimipisteessä on
% oltava. Samaa toleranssia käytetään myös kun tarkastetaan tuleeko
% sakkokierros.
% 
% 
% Gradientti ja hessen matriisi määritetään omien funktioidensa avulla
% numeerisilla aproksimaatioilla. Tiedoston lopusta löytyy funktiot numGrad
% ja hesse. 
% 
% Sakkofunkion muodostamisessa hyödynnetään max funktiota, joka palauttaa
% arvoista 0 ja g(xk) aina suuremman. Täten, kun ehdot ovat muotoa 
% g(xk) >=0, niin jos rajoite toteutuu saadaan sakoksi 0 (arvo silloin
% negatiivinen) ja jos ei niin g(xk)^2.



function [xk,ite,sakkokierros] = ngs(f,g,xk,tol)

    % Arvojen alustus 
    s = 1;
    alpha = 0.5;
    beta = 0.5;
    ite = 0;
    sakkokierros  = 0;
    sakko = @(x) sum(max(0,g(x)).^2);
    P = @(x,s) f(x) + s*sakko(x);
    Pg = numGrad(P,xk,s);
    H = hesse(P,xk,s);

    while norm(Pg) > tol

        ite = ite + 1;
        % Lopetetaan, jos tulee iteraatioita täyteen tietyn verran.
        if ite > 1e+3
            break
        end

        % Suunnan valinta
        om = eig(H);
        if all(om > 0)
            d = -H\Pg;
        else
            d = -Pg;
        end

        % Viivahaku
        a=1;
        while P(xk+a*d,s) > P(xk,s)+alpha*a*Pg'*d
            a = beta*a;
        end
        
        % Arvojen päivitys
        xk = xk+ a*d;
        Pg = numGrad(P,xk,s);
        H = hesse(P,xk,s);

        % Rajoitteiden toteutumisen tarkastus
        if sakko(xk) > tol
            s = s*2;
            sakkokierros = sakkokierros +1;
        end  
    end
end



function g = numGrad(f,x,s)
    h = 1e-3;
    n = length(x);
    g = zeros(n,1);
    for i = 1:n
        xi1 = x;
        xi2 = x;
        xi1(i) = x(i)+h;
        xi2(i) = x(i)-h;
        g(i) = (f(xi1,s)-f(xi2,s))/(2*h);
    end
end


function H = hesse(f,x,s)
    h = 1e-3;
    n = length(x);
    H = zeros(n);
    d= h*eye(n);
    for i = 1:n
        for j = 1:n
            H(i,j) = (...
                f(x+d(:,i)+d(:,j),s)...
                -f(x+d(:,i)-d(:,j),s)...
                -f(x-d(:,i)+d(:,j),s)...
                +f(x-d(:,i)-d(:,j),s))/(4*h^2);
        end
    end   
end




