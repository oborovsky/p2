% разностная схема для уравнения
% -m*D2(u) + B(x)*D1(u) + C(x)*u = F(x)
% u(0) = 0 , u(1) = 1; при малых m и С(x)>= 0; B(x) < 0;
% рассмотрена аппроксимация уравнения
% -m*D2(u) - D1(u) = 0;
% с решением u(x)  = (exp(-x/m)-1)/(exp(-1/m)-1);

mm = 0.5:-0.01:0.01;
for j = 1:length(mm)
    N = 50;
    for i=1:2

        N = 2*N;
        h = 1/N;
        m = mm(j)
        a = 0.15;
        s = 0:h:1;
        x1 = X1(s,m,a);
        x2 = X2(s,m,a);
    %     x = s;
    %     subplot(1,2,1);
    %     plot(s,x1,'r');
    %     hold on;
    %     plot(s,x,'b');
    %     hold off;

        [m1,d1] = makeLinearSystem2(x1,m);
        [m2,d2] = makeLinearSystem3(x2,m);
        u1(i) = {shufle(m1,d1)};
        y1(i) = {U(x1,m)};
        u2(i) = {shufle(m2,d2)};
        y2(i) = {U(x2,m)};
        ee1 = u1{i} -y1{i};
        ee2 = u2{i} - y2{i};
        e1(i) = max(abs(ee1));
        e2(i) = max(abs(ee2));
        if  i >= 2
            fprintf('e1(%d)=%1.15f,e1(%d)=%1.15f, p1 = %1.15f\n',i-1,e1(i-1),i,e1(i), abs(log2(e1(i-1)/e1(i))));
            fprintf('e2(%d)=%1.15f,e2(%d)=%1.15f, p2 = %1.15f\n',i-1,e2(i-1),i,e2(i), abs(log2(e2(i-1)/e2(i))));
        end
    %     subplot(1,2,2);
%         figure(1);
%         plot(x1,y1{i},'r',x1,u1{i},'g');
%         figure(2);
%         plot(x2,y2{i},'r',x2,u2{i},'g');
    %     plot(x,u{i},'g');

    end
    eps1(j) = e1(length(e1));
    eps2(j) = e2(length(e2));
end
 semilogy(mm,eps1,'r');
 title('зависимость погрешности от m');
%  hold on;
%  semilogy(mm,eps2,'g');
%  hold off;
 

function [m,d] = makeLinearSystem2(x,m)
% аппроксимация 2 - го порядка
     Hm = left(x); 
     Hp = right(x);
     
    aa = [0];
    bb = [1];
    cc = [0];
    dd = [0];
    
    l = length(x);
    
    for i=2:l-1
        xi = x(i);
        hm = Hm(i);
        hp = Hp(i);
        s = hm + hp;
        d = hp - hm;
        p = hm*hp;
        
        bp = B(xi + hp/2);
        bm = B(xi - hm/2);
        b = B(xi);
        cm = C(xi-hm);
        cp = C(xi+hp);
        c = C(xi);
        
        k1 = -2*m/(s*hm) - b*hp/(s*hm) + 2*d*bm/(3*s*hm) - d*hp*cm/(3*s*hm);
        k3 = -2*m/(s*hp) + b*hm/(s*hp) + 2*d*bp/(3*s*hp) + d*hm*cp/(3*s*hp);
        k2 = 2*m/p + b*d/p - 2*d*(bp*hm + bm*hp)/(3*s*p) + (1+d/p)*c;
        
        aa(i) = k1;
        bb(i) = k2;
        cc(i) = k3;
        
        b3 = d*hm/(3*s*hp);
        b1 = -d*hp/(3*s*hm);
        b2 = (1+d/p);
        dd(i) = b1*F(xi - hm) + b2*F(xi) + b3*F(xi+hp);
    end
    aa(l) = 0;
    bb(l) = 1;
    cc(l) = 0;
    dd(l) = 1;
    m = [aa;bb;cc];
    d = dd';
end
function [m,d] = makeLinearSystem3(x,m)
% аппроксимация 3-го порядка
    Hm = left(x); 
     Hp = right(x);
     
    aa = [0];
    bb = [1];
    cc = [0];
    dd = [0];
    
    l = length(x);
    
    for i=2:l-1
        xi = x(i);
        hm = Hm(i);
        hp = Hp(i);
        s = hm + hp;
        d = hp - hm;
        p = hm*hp;
        
        bp = B(xi + hp/2);
        bm = B(xi - hm/2);
        b = B(xi);
        db = -hp*B(xi-hm)/(s*hm) + d*B(xi)/p + hm*B(xi+hp)/(s*hp);
        lb = 2*B(xi-hm)/(s*hm) - 2*B(xi)/p + 2*B(xi+hp)/(s*hp);
        
        cm = C(xi-hm);
        cp = C(xi+hp);
        c = C(xi);
        
        g1 = -m + p*db/6;
        g2 = b + p*lb/12;
        r = d/3 - b*(p-d^2/3)/(12*m);
        
        k1 = 2*g1/(s*hm) - g2*hp/(s*hm) - 2*r*bm/(s*hm) + (-r*hp/(s*hm) + (p+d^2)/(6*s*hm))*cm;
        k3 = 2*g1/(s*hp) + g2*hm/(s*hp) - 2*r*bp/(s*hp) + ( r*hm/(s*hp) + (p+d^2)/(6*s*hp))*cp;
        k2 = -2*g1/p + g2*d/p + 2*r*(bp*hm + bm*hp)/(s*p) + (1 + d*r/p - (p+d^2)/(6*p))*c;
        
        aa(i) = k1;
        bb(i) = k2;
        cc(i) = k3;
        
        b3 = r*hm/(s*hp) + (p+d^2)/(6*s*hp);
        b1 = -r*hp/(s*hm) + (p+d^2)/(6*s*hm);
        b2 = 1 + d*r/p - (p+d^2)/(6*p);
        dd(i) = b1*F(xi - hm) + b2*F(xi) + b3*F(xi+hp);
    end
    aa(l) = 0;
    bb(l) = 1;
    cc(l) = 0;
    dd(l) = 1;
    m = [aa;bb;cc];
    d = dd';
end
function y=B(x)
% коэфф. при D1u
    y  = -1; %(1)(2)
end
function y=C(x)
% коэфф. при неизвестной функции
    y = 0; %(1)(2)
end
function y=F(x)
% правая часть уравнения
    y = 0; %(1)
%     y = x.^2; %(2)
end
function y=U(x,m)
% решение уравнения
    y = (exp(-x/m) - 1)/(exp(-1/m) - 1); %(1)
%       y = (exp(-x/m) - exp(-1/m))/(1- exp(-1/m));
end
function y=X1(s,m,a)
% X1 - сетка для схемы 2-го прорядка
    p = 2/(1+m^a);
    l = length(s);
    
    for i = 1:l
        x = s(i);
        if x >= 0 & x <= 1/2
            y(i) = m*(1/(1-p*x) - 1);
        else
            K1 = 4*(1-m^(1-a)-p*m^(1-2*a)*(1+m^a)/2);
            y(i) = m^(1-a) + p*m^(1-2*a)*(1+m^a)*(x-1/2) + K1*(x-1/2)^2;
        end
    end
end
function y=X2(s,m,a)
%X2 - сетка (с гладкостью 2-го порядка) для схемы 3-го порядка
    p = 2/(1+m^a);
    l = length(s);
    
    for i = 1:l
        x = s(i);
        if x >= 0 & x <= 1/2
            y(i) = m*(1/(1-p*x) - 1);
        else
            K1 = p^2*m/(1-p/2)^3;
            K2 = 8 - 8*m^(1-a) - 4*p*m^(1-2*a)*(1+m^a) - 2*K1;
            y(i) = m^(1-a) + p*m^(1-2*a)*(1+m^a)*(x-1/2) + K1*(x-1/2)^2 + K2*(x-1/2)^3;
        end
    end
end
function y = left(x)
% left - вычисляет вектор h-;
    l = length(x);
    y(1) = 0;
    for i = 2:l
        y(i) = x(i) - x(i-1);
    end
end
function y = right(x)
% right - вычисляет вектор h+;
    l = length(x);
    for i = 1:l-1
        y(i) = x(i+1) - x(i);
    end
    y(l) = 0;
end
function x = shufle(m,d)
% SHUFLE реализация метода прогонки для уравнения:
% a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i), i = 1,...,n
% a(1)=c(n)= 0;
% аргументы :
% m - массив (3хn), где
% m(1,:) = a(1:end); m(2,:) = b(1:end); m(3,:) = c(1:end)
% d - столбец правых частей
    a = m(1,:);
    b = m(2,:);
    c = m(3,:);
    
    s = length(b);
    
    et(1) = 0;
    ks(1) = 0;
    
    for i=1:s
        ks(i+1) = (-c(i)) / (a(i)*ks(i) + b(i));
        et(i+1) = (d(i) - a(i) * et(i)) / (a(i) * ks(i) + b(i));
    end
    
    x(s+1) = 0;
    
    for j=0:s-1
        x(s-j) = ks(s-j+1)*x(s-j+1) + et(s-j+1);
    end
    x(s+1) = [];
end