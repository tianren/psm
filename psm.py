#!/usr/bin/python3
__author__ = 'Tianren Liu'

import sys, os
from math import gcd
import numpy as np
# import sympy
# try:
#     import sympy
# except ImportError:
#     sympy = None

# modular multiplicative inverse
# mod: assume prime
def modpow(a,k,mod):
    res = 1
    curbit = a
    while k > 0:
        if k%2 == 1:
            res = (res * curbit) % mod
        curbit = (curbit * curbit) % mod
        k //= 2
    return res

# modular multiplicative inverse
# mod: assume prime
def modinv(a,mod):
    return modpow(a,mod-2,mod) % mod

# check if matrix M span row vector (1,0,...,0)
def span_eigen(M,mod):
    M = M.copy()
    M %= mod
    h,w = M.shape
    C = np.identity(h,dtype=int)

    lead = 0
    for j in range(w):
        # print (j,lead,":::")
        # sys.stdin.readline()
        # # print (M)
        # # search non-zero
        for k in range(lead, h):
            if M[k,j] != 0:
                if k != lead:
                    C[[lead, k]] = C[[k, lead]]
                    M[[lead, k]] = M[[k, lead]]
                C[lead] *= modinv(M[lead,j], mod)
                C[lead] %= mod
                M[lead] *= modinv(M[lead,j], mod)
                M[lead] %= mod
                for i in range(h):
                    if i != lead:
                        C[i] -= M[i,j] * C[lead]
                        C[i] %= mod
                        M[i] -= M[i,j] * M[lead]
                        M[i] %= mod
                lead += 1
                break
        else:
            if j==0 or M[0,j] != 0:
                return None
    return C[0]

def factorial(k):
    res = 1
    for _k in range(1,k+1):
        res *= _k
    return res


# enumerate sub-squences,
# par: a squence, assuming non-increasing
# maxsum: upper bound on sub-squence sum
def enum_subsquence(par):
    if len(par) == 0:
        yield tuple()
    else:
        tail = par[-1]
        tailstart = par.index(tail)
        tailcount = len(par) - tailstart
        for _head in enum_subsquence(par[:tailstart]):
            yield _head
            for i in range(tailcount):
                _head = _head + (tail,)
                yield _head
def enum_subsquence_with_complement(par):
    if len(par) == 0:
        yield tuple(), tuple()
    else:
        tail = par[-1]
        tailstart = par.index(tail)
        tailcount = len(par) - tailstart
        for _head,_comp in enum_subsquence_with_complement(par[:tailstart]):
            for i in range(tailcount+1):
                yield _head + i*(tail,), _comp + (tailcount-i)*(tail,)


class partition:
    def __init__(self, n, maxs=None, mins=1):
        if maxs == None:
            maxs = n
        self.space = n
        self.shape_head, self.shape_tail = maxs, mins
    def nextshape(self,x):
        return x-1
    def fit(self, shape, space):
        return shape <= space
    def subtract(self, shape, space):
        return space - shape
    def add(self, shape, space):
        return space + shape
    def empty(self, space):
        return space == 0
class partition_twosets:
    def __init__(self, n, maxs=None):
        if maxs == None:
            maxs = 2*n
        self.space,self.n,self.maxs = (n,n), n, maxs
        self.shape_head, self.shape_tail = (min(n,maxs), maxs-min(n,maxs)), (0,1)
    def nextshape(self,x):
        if x[1] > 0:
            return (x[0], x[1]-1)
        else:
            return (x[0]-1, min(self.n, self.maxs-(x[0]-1)))
    def fit(self, shape, space):
        return shape[0] <= space[0] and shape[1] <= space[1]
    def subtract(self, shape, space):
        return (space[0] - shape[0], space[1] - shape[1])
    def add(self, shape, space):
        return (space[0] + shape[0], space[1] + shape[1])
    def empty(self, space):
        return space == (0,0)

def pattern_iter(P):
    fixed, freespace, shape = list(), P.space, P.shape_head
    while True:
        # print (fixed, freespace, shape)
        if P.fit(shape, freespace):
            fixed.append(shape)
            freespace = P.subtract(shape, freespace)
        else:
            if P.empty(freespace):
                yield tuple(fixed)
                shape = fixed.pop()
                freespace = P.add(shape, freespace)
            while shape == P.shape_tail and len(fixed) > 0:
                shape = fixed.pop()
                freespace = P.add(shape, freespace)
            if shape != P.shape_tail:
                shape = P.nextshape(shape)
            else:
                break

class psm_2blocks_per_party(partition):
    def __init__(self, n, CC = None, _CC = lambda x:x-1, easytermmax=None, _easytermmax=lambda x:x+1):
        if CC is None:
            CC = _CC(n)
        if easytermmax is None:
            easytermmax = _easytermmax(n)

        self.easytermmax = easytermmax
        partition.__init__(self, 2*n, CC)
        self.target = tuple()
    def combination_num(self, shapes):
        n = sum(shapes)
        res = factorial(n)
        p,c = 0,0
        for s in shapes:
            res //= factorial(s)
            if s == p:
                c += 1
                res //= c
            else:
                p,c = s,1
        return res
    def iseasyterm(self, shapes):
        return sum(shapes) <= self.easytermmax
    def str_paddedterm(self, shape):
        return "Σf(barX{})".format(shape)
    def str_hardterm(self, shape):
        return "Σf(R{})".format(shape)
    # def str_paddedterm(self, shape):
    #     return "\sumFbarX{{{}}}".format(','.join((str(_) for _ in shape)))
    # def str_hardterm(self, shape):
    #     return "\sumFR{{{}}}".format(','.join((str(_) for _ in shape)))

class psm_2party_tradeoff(partition):
    def __init__(self, K, CC, otherCC=None):
        if otherCC is None:
            otherCC = K - CC
        self.otherCC = otherCC
        partition.__init__(self, K, CC)
        self.target = tuple()
    combination_num = psm_2blocks_per_party.combination_num
    def iseasyterm(self, shapes):
        return sum(shapes) <= self.otherCC
    def str_paddedterm(self, shape):
        return "Σf(barX{},Z)".format(shape)
    def str_hardterm(self, shape):
        return "Σf(R{},Z)".format(shape)

class psm_2party(partition_twosets):
    def __init__(self, K, CC):
        self.CC = CC
        partition_twosets.__init__(self, K, CC)
        self.target = tuple()
    def combination_num(self, shapes):
        res = factorial(sum((s[0] for s in shapes))) * factorial(sum((s[1] for s in shapes)))
        p,c = None, 0
        for s in shapes:
            res //= factorial(s[0]) * factorial(s[1])
            if s == p:
                c += 1
                res //= c
            else:
                p,c = s,1
        return res
    def iseasyterm(self, shapes):
        return sum((s[0] for s in shapes)) <= self.CC or sum((s[1] for s in shapes)) <= self.CC
    def str_paddedterm(self, shape):
        return "Σf(barX{})".format(shape)
    def str_hardterm(self, shape):
        return "Σf(R{})".format(shape)

def general_solver(P, mod, verbose=False):
    _mod = (lambda x:x) if mod is None else (lambda x:x%mod)

    rterms = tuple(pattern_iter(P))
    freeterms = list()
    coeffmatrix = list()

    for rterm in rterms:
        coeffvec = len(freeterms) * [0,]
        for s,c in enum_subsquence_with_complement(rterm):
            if not P.iseasyterm(c):
                try:
                    si = freeterms.index(s)
                    coeffvec[si] = _mod(P.combination_num(c))
                except ValueError:
                    freeterms.append(s)
                    coeffvec.append(_mod(P.combination_num(c)))
        coeffmatrix.append(coeffvec)
        if verbose:
            print ("{} = easy terms + {}".format(P.str_paddedterm(rterm), " + ".join(("{}*{}".format(coeff,P.str_hardterm(freeterms[si])) for si,coeff in enumerate(coeffvec) if coeff != 0))))

    M = np.zeros((len(rterms), len(freeterms)), dtype=int)
    for _i, vec in enumerate(coeffmatrix):
        M[_i,:len(vec)] = vec
    if verbose:
        print(M)
        print(M.shape)

    algo = span_eigen(M,mod)
    if algo is not None:
        print ("{} = easy terms".format(P.str_hardterm(tuple())) + "".join((" + {}*{}".format(t,P.str_paddedterm(rterms[i])) for i,t in enumerate(algo) if t != 0)) + " mod {}".format(mod))
    return algo

if __name__ == "__main__":
    try:
        np.set_printoptions(linewidth=os.get_terminal_size().columns)
    except OSError:
        pass

    # Search for a k-party PSM protocol
    general_solver(psm_2blocks_per_party(7), mod=100000000003, verbose=True)
    # general_solver(psm_2blocks_per_party(7), mod=19, verbose=True)
    # general_solver(psm_2blocks_per_party(8), mod=19, verbose=True)
    # general_solver(psm_2blocks_per_party(9), mod=19, verbose=True)
    # general_solver(psm_2blocks_per_party(10), mod=11, verbose=True)
    # general_solver(psm_2blocks_per_party(11), mod=23, verbose=True)
    # general_solver(psm_2blocks_per_party(11), mod=11, verbose=True)
    # general_solver(psm_2blocks_per_party(12), mod=13, verbose=True)
    # general_solver(psm_2blocks_per_party(13), mod=29, verbose=True)
    # general_solver(psm_2blocks_per_party(14), mod=29, verbose=True)
    # general_solver(psm_2blocks_per_party(15), mod=2, verbose=True)
    # general_solver(psm_2blocks_per_party(16), mod=17, verbose=True)
    # general_solver(psm_2blocks_per_party(17), mod=71, verbose=True)
    # general_solver(psm_2blocks_per_party(17), mod=17, verbose=True)
    # general_solver(psm_2blocks_per_party(18), mod=19, verbose=True)
    # general_solver(psm_2blocks_per_party(19), mod=71, verbose=True)
    # general_solver(psm_2blocks_per_party(20), mod=3001, verbose=True)

    # Search for a 2-party PSM with unbalanced c.c.
    # for i in range(1,12):
    #    print (i)
    #    general_solver(psm_2party_tradeoff(12,i), mod=13)
    # for i in range(1,20):
    #    print (i)
    #    general_solver(psm_2party_tradeoff(20,i), mod=71)
