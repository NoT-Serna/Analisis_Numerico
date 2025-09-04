function s = Simpson(a,b,n,f)
    h = (b-a)/n;
    xl = f(a)+f(b);
    xl1 = 0;
    xl2= 0;
    for i = 1:n-1
        x = a+i*h;
        if mod(i,2)==0
            xl2 = xl2 + f(x);
        else
            xl1 = xl1 + f(x);
        end
    end
    s = (h/3) * (xl + 4*xl1 + 2*xl2);
end
