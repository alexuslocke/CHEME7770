syms u v a n

dydt1 = -u + a/(1+v^n)
dydt2 = -v + a/(1+u^n)

jacobian([dydt1, dydt2], [u v])
