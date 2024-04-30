import rtcore
from numpy import *


ps = array([[1., 2.2, 8.5], [-5.1, 0.7, 3.5], [-1., 0., 7.]])
dr = array([[7.9, 6.8, 4.2], [4.1, 7.7, 5.0], [4., -3., 6.]])

freq = 1e7
dcr = 45.62362
ris = 25.0

res = empty(2, dtype=double)
tb =  empty(2, dtype=double)
od =  empty(2, dtype=double)
flg = empty(2, dtype=short)
nw = 1;

print 'Entering trace_beam() ...'

rtcore.trace_beam(ps, dr, freq, dcr, ris, res, tb, od, flg, nw)

