syms x R c0 c2 c4
f(x) = (sqrt(1-(x/R)^2))*(c0+c2*(x/R)^2+c4*(x/R)^4);
pretty(f(x))
df(x) = diff(f,x);
pretty(simplify(df(x)))

R=3.91;
c0=0.81;
b=4.844/2;
h=2.58;
eq1=vpa(subs(f(b))==h)
eq2=vpa(subs(df(b))==0)

[A,B] = equationsToMatrix([eq1, eq2], [c2, c4])
X = linsolve(A,B)