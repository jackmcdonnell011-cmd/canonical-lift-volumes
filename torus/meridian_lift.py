"""Exact meridians of the lifted knots, by lifting mu_K through the covering.

SnapPy's cover() numbers the tetrahedra of a degree-d cover as  i = t*d + s
(base tetrahedron t, sheet s) with the same vertex labelling and gluing
permutations as the base (checked below).  The full preimage of the base
meridian mu_K is therefore obtained by copying, into every sheet, the base
meridian's normal-curve data at the vertices lying on the K cusp.

Because rho(mu_K) = 1, every lift of mu_K is closed, so on a K-lift cusp of
degree m the preimage is m parallel copies of that cusp's meridian.  We write
the preimage into the cover's triangulation file in place of SnapPy's meridian
(keeping its longitude), reload, and compare cusp shapes: holonomy is linear in
the curve data, so
    H(preimage) / H(mu_s) = tau_old / tau_new = m * (p + q*tau_old)
where (mu_s, lambda_s) is SnapPy's basis on that cusp.  (p, q) must come out as
integers (to 1e-6) and coprime; that is the exact check that the slope we fill
is the lift of mu_K.
"""
import snappy


def _parse(text):
    lines = text.splitlines()
    last = max(k for k, l in enumerate(lines) if 'torus' in l or 'Klein' in l)
    i = last + 1
    while not lines[i].strip():
        i += 1
    ntet = int(lines[i].strip())
    head = lines[:i + 1]
    tets = []
    j = i + 1
    for _ in range(ntet):
        nb = lines[j]; pm = lines[j + 1]; cu = lines[j + 2]
        curves = [list(map(int, lines[j + 3 + r].split())) for r in range(4)]
        sh = lines[j + 7]
        tets.append(dict(nb=nb, pm=pm, cu=list(map(int, cu.split())), curves=curves, sh=sh))
        j += 8
        while j < len(lines) and not lines[j].strip():
            j += 1
    return head, tets


def _unparse(head, tets):
    out = list(head)
    for t in tets:
        out += [t['nb'], t['pm'], ' ' + ' '.join(f'{c:4d}' for c in t['cu'])]
        for r in t['curves']:
            out.append(' ' + ' '.join(f'{x:2d}' for x in r))
        out.append(t['sh'])
        out.append('')
    return '\n'.join(out) + '\n'


def lifted_meridian_slopes(M, cuspK, C, d, Kcusps_deg, hint=None):
    """M: base manifold, C: its degree-d cover. Kcusps_deg: {cover cusp: degree m}.
    Returns {cover cusp: (p, q)} and a dict of diagnostics."""
    bh, bt = _parse(M._to_string())
    ch, ct = _parse(C._to_string())
    assert len(ct) == len(bt) * d
    # check the tetrahedron correspondence i = t*d + s
    for i, T in enumerate(ct):
        B = bt[i // d]
        assert T['pm'].split() == B['pm'].split(), 'gluing permutations differ'
        bnb = list(map(int, B['nb'].split()))
        cnb = list(map(int, T['nb'].split()))
        assert all(cn // d == bn for cn, bn in zip(cnb, bnb)), 'neighbour tetrahedra differ'
    # Put the lifted curve in the meridian slot, unless the expected answer is SnapPy's
    # current longitude (0,1): then use the longitude slot instead, so that the file never
    # contains two parallel peripheral curves (that crashes the SnapPea kernel).
    hint = hint or {}
    slot = {c: (2 if hint.get(c) == (0, 1) else 0) for c in Kcusps_deg}
    for i, T in enumerate(ct):
        B = bt[i // d]
        for v in range(4):
            if B['cu'][v] == cuspK:
                c = T['cu'][v]
                for r in (0, 1):  # right- and left-handed sheets of the base meridian
                    T['curves'][slot[c] + r][4 * v:4 * v + 4] = B['curves'][r][4 * v:4 * v + 4]
    text = _unparse(ch, ct)
    N = snappy.Manifold(text)
    old = [complex(c['shape']) for c in C.cusp_info()]
    new = [complex(c['shape']) for c in N.cusp_info()]
    slopes, diag = {}, {}
    for c, m in Kcusps_deg.items():
        tau, tau2 = old[c], new[c]
        from math import gcd
        best = None
        # SnapPea computes a curve's holonomy by tracing it, so for a degree-m cusp the
        # value may be that of ONE lifted component (divisor 1) or of all m (divisor m).
        for div in (1, m):
            # meridian slot: H(lift)/H(mu_s) = tau/tau2 ; longitude slot: H(lift)/H(mu_s) = tau2
            z = ((tau / tau2) if slot[c] == 0 else tau2) / div
            q = z.imag / tau.imag
            p = z.real - q * tau.real
            P, Q = round(p), round(q)
            err = abs(p - P) + abs(q - Q)
            if err < 1e-6 and gcd(abs(P), abs(Q)) == 1:
                best = (P, Q, err, div)
                break
        if best is None:
            slopes[c] = None
            diag[c] = dict(ok=False)
            continue
        P, Q, err, div = best
        if P < 0 or (P == 0 and Q < 0):
            P, Q = -P, -Q
        slopes[c] = (P, Q)
        diag[c] = dict(err=err, ok=True, divisor=div)
    return slopes, diag
