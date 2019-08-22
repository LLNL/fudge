x = ( a + b ) / 2
d = ( b - a ) / 2.


f[y_] = Erf[y] / 2.
f[y_] = - Exp[-y^2] / 2. / Sqrt[Pi]
f[y_] = 0.25 * Erf[y] - y Exp[-y^2] / 2. / Sqrt[Pi]
f[y_] = -( y^2 + 1 ) Exp[-y^2] / 2. / Sqrt[Pi]
f[y_] = 3. Erf[y] / 8. - y ( 2. y^2 + 3. ) Exp[-y^2] / 4. / Sqrt[Pi]

Exp[ -x^2 - d^2 ] Sinh[ 2. x d]

f[b] - f[a]
