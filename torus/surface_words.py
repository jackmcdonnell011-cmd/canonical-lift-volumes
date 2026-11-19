"""Words in pi_1(S) = Gamma, for S = H/Gamma with Gamma torsion-free of finite index in PSL2(Z).

Free basis: the Farey triangulation of S is dual to a trivalent ribbon graph (vertices = triangles =
orbits of R = S.X on the cosets, edges = orbits of S).  Take a spanning tree of that graph; each edge
{x, x.S} not in the tree gives one generator P(x) S P(x.S)^-1, where P(y) is the tree path from coset 0
to coset y.  These are the side pairings of the ideal polygon made of the triangles, and there are
d/2 - d/3 + 1 = d/6 + 1 = 2g + n - 1 of them.  Generator i is written i+1, its inverse -(i+1).

A closed geodesic on S is a conjugacy class of Gamma, i.e. a cyclic word in this basis up to rotation.
"""
from collections import deque
import psl2

R = psl2.mul(psl2.S, psl2.X)


class Basis:
    def __init__(self, G):
        self.G = G
        d = G.d
        pS, pR = G.perm(psl2.S), G.perm(R)
        self.pS, self.pR = pS, pR
        # vertices of the ribbon graph = R-orbits; P[y] = matrix taking coset 0 to coset y along the tree
        P = {0: psl2.I2}
        for y, M in ((pR[0], R), (pR[pR[0]], psl2.mul(R, R))):
            P[y] = M
        tree = set()
        q = deque([0])
        while q:
            x0 = q.popleft()
            for x in (x0, pR[x0], pR[pR[x0]]):      # the three darts of this vertex
                y = pS[x]
                if y in P:
                    continue
                tree.add(min(x, y))
                P[y] = psl2.mul(P[x], psl2.S)
                P[pR[y]] = psl2.mul(P[y], R)
                P[pR[pR[y]]] = psl2.mul(P[y], psl2.mul(R, R))
                q.append(y)
        assert len(P) == d
        self.P = P
        self.edge_gen = {}
        self.gens = []
        for x in range(d):
            y = pS[x]
            if x < y and x not in tree:
                g = psl2.mul(psl2.mul(P[x], psl2.S), psl2.inv(P[y]))
                assert G.act(0, g) == 0
                self.edge_gen[x] = len(self.gens) + 1      # crossing x -> x.S is +gen, y -> x is -gen
                self.gens.append(g)
        assert len(self.gens) == d // 6 + 1, (len(self.gens), d)

    def _cross(self, x, out):
        """cross the edge at dart x; record a generator if it is not a tree edge."""
        y = self.pS[x]
        e = min(x, y)
        if e in self.edge_gen:
            out.append(self.edge_gen[e] if x == e else -self.edge_gen[e])
        return y

    def word(self, M):
        """Free-reduced word in the basis for M in Gamma."""
        assert self.G.act(0, M) == 0
        pR = self.pR
        x, out = 0, []
        for t, k in psl2.st_word(M):
            if t == 'S':
                x = self._cross(x, out)
            elif k > 0:
                for _ in range(k):                       # X = S.R
                    x = pR[self._cross(x, out)]
            else:
                for _ in range(-k):                      # X^-1 = R^2.S
                    x = self._cross(pR[pR[x]], out)
        assert x == 0
        return reduce(out)

    def matrix(self, word):
        M = psl2.I2
        for i in word:
            g = self.gens[abs(i) - 1]
            M = psl2.mul(M, g if i > 0 else psl2.inv(g))
        return M


def reduce(w):
    out = []
    for i in w:
        if out and out[-1] == -i:
            out.pop()
        else:
            out.append(i)
    return out


def cyclic_reduce(w):
    w = reduce(w)
    while len(w) > 1 and w[0] == -w[-1]:
        w = w[1:-1]
    return w


def cyclic_canon(w):
    """canonical rotation of a cyclically reduced word (an oriented closed curve)."""
    w = cyclic_reduce(w)
    n = len(w)
    if n == 0:
        return ()
    return min(tuple(w[i:] + w[:i]) for i in range(n))


def inverse(w):
    return [-i for i in reversed(w)]


def unoriented_canon(w):
    return min(cyclic_canon(w), cyclic_canon(inverse(w)))


def is_primitive_cyclic(w):
    n = len(w)
    return all(not (n % p == 0 and w[:p] * (n // p) == list(w)) for p in range(1, n))


def letters(w, names='abcdefghijklmnop'):
    return ''.join(names[i - 1] if i > 0 else names[-i - 1].upper() for i in w)


def enumerate_curves(r, N):
    """One cyclically reduced primitive word per unoriented closed curve, word length <= N."""
    seen = set()
    out = []
    def rec(w):
        if w:
            if w[0] != -w[-1] and is_primitive_cyclic(w):
                c = unoriented_canon(w)
                if c not in seen:
                    seen.add(c)
                    out.append(list(c))
        if len(w) == N:
            return
        for i in list(range(1, r + 1)) + list(range(-r, 0)):
            if w and w[-1] == -i:
                continue
            rec(w + [i])
    rec([])
    return out
