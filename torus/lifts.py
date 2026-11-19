"""Canonical-lift complements on S = H^2/Gamma via finite covers of the modular link.

M_w = S^3 - (trefoil u K_w) = UT(Mod) - (canonical lift of w)  (get_dt from modular_dt.py)
For Gamma <= PSL2(Z) torsion-free of index d, the cover of M_w given by
    pi_1(M_w) --rho--> PSL2(Z) --> Sym(Gamma\\PSL2(Z))
is UT(H/Gamma) minus the full preimage of the lift of w, i.e. minus the canonical
lifts of all closed geodesics on S lying over w (one per cycle of w on the cosets).
Filling a K-lift cusp along its meridian puts that knot back.
"""
import math
import snappy
import psl2
import rho as rhomod
from modular_dt import get_dt


def word_length(w):
    t = abs(psl2.tr(psl2.xy_matrix(w)))
    return 2 * math.acosh(t / 2)


def _inv(p):
    q = [0] * len(p)
    for i, j in enumerate(p):
        q[j] = i
    return q


def _act_word(perms, word, x):
    for ch in word:
        p = perms[ch] if ch.islower() else _inv(perms[ch.lower()])
        x = p[x]
    return x


def _cycles(p):
    seen, out = set(), []
    for i in range(len(p)):
        if i in seen:
            continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j)
            cyc.append(j)
            j = p[j]
        out.append(cyc)
    return out


def solve_rho(dt, cuspK=0, cuspT=1, tries=80, attempts=None, cont_tries=None, hard=None):
    import os as _os
    if hard is None:
        hard = _os.environ.get('GN_HARD', '') == '1'
    attempts = attempts or (12 if hard else 6)
    cont_tries = cont_tries or (30 if hard else 15)
    """Robust solve on fresh (re-loaded, hence reproducible) triangulations:
    fast LM first, then complex-detour continuation (rho_cont)."""
    import rho_cont
    import rho_exact
    # exact: fill K, send the trefoil group <a, b | a^2 b^3> to (S, R), carry back (see rho_exact.py)
    M = snappy.Manifold(snappy.Link('DT: ' + dt).exterior()._to_string())
    gm = rho_exact.find_rho_exact(M, cuspK, cuspT)
    if gm is not None:
        return gm, M, M.fundamental_group(), 'exact', 0
    for attempt in range(attempts):
        M = snappy.Manifold(snappy.Link('DT: ' + dt).exterior()._to_string())
        G = M.fundamental_group()
        gm = rhomod.find_rho(G, cuspK, cuspT, tries=tries, seed=attempt)
        if gm is not None:
            return gm, M, G, 'LM', attempt
        gm = rho_cont.find_rho_cont(G, cuspK, cuspT, theta0=math.pi / 3, tries=cont_tries, seed=attempt)
        if gm is not None:
            return gm, M, G, 'continuation', attempt
        if hard:
            import rho_fast
            gm = rho_fast.find_rho_fast(G, cuspK, cuspT, seed=attempt, direct_tries=20, cont_tries=40)
            if gm is not None:
                return gm, M, G, 'newton-continuation', attempt
    return None, None, None, None, attempts


class ModularLink:
    """The modular link complement for word w, plus the Ghys map rho."""

    def __init__(self, w, tries=80, cache_dir=None):
        self.w = w
        self.dt = get_dt(w)
        if self.dt is None:
            raise ValueError('word not primitive/canonical: ' + w)
        self.cuspK, self.cuspT = 0, 1
        gm = None
        # Link.exterior() randomises the triangulation (and hence the presentation);
        # the numerical solve is easier for some presentations, so retry on fresh ones.
        self.solver = None
        gm = None
        cache = None
        if cache_dir:
            import json, os
            fn = os.path.join(cache_dir, w + '.json')
            if os.path.exists(fn):
                cache = json.load(open(fn))
        attempt = 0
        if cache is not None:
            self.M = snappy.Manifold(cache['tri'])
            self.G = self.M.fundamental_group()
            gm0 = {g: tuple(v) for g, v in cache['gm'].items()}
            if rhomod.verify(self.G, gm0, self.cuspK, self.cuspT):
                gm = gm0
                self.solver = 'cache:' + cache.get('solver', '?')
        if gm is None:
            gm, self.M, self.G, self.solver, attempt = solve_rho(self.dt, self.cuspK, self.cuspT, tries)
            if cache_dir and gm is not None:
                import json, os
                os.makedirs(cache_dir, exist_ok=True)
                fn = os.path.join(cache_dir, w + '.json')
                tmp = fn + '.tmp'
                json.dump(dict(word=w, dt=self.dt, tri=self.M._to_string(), gm=gm, solver=self.solver,
                               relators=self.G.relators()), open(tmp, 'w'))
                os.replace(tmp, fn)
        if gm is None:
            raise RuntimeError('rho not found for ' + w)
        self.attempts = attempt + 1
        per = self.G.peripheral_curves()
        lK = psl2.word_matrix(per[self.cuspK][1], gm)
        target = psl2.canon_word(w)
        cl = {psl2.canon_class(lK), psl2.canon_class(psl2.inv(lK))}
        clJ = {psl2.canon_class(psl2.conjJ(lK)), psl2.canon_class(psl2.inv(psl2.conjJ(lK)))}
        if target in cl:
            self.twisted = False
            self.gm = gm
        elif target in clJ:
            self.twisted = True
            self.gm = rhomod.twist(gm)
        else:
            raise RuntimeError(f'lambda_K class {cl} does not match word {w}')
        lK = psl2.word_matrix(per[self.cuspK][1], self.gm)
        # orientation of K relative to the word: +1 if lambda_K ~ w, -1 if ~ w^{-1}
        self.Korient = 1 if psl2.canon_class(lK) == target else -1
        self.lK = lK
        self.vol = float(self.M.volume())
        self.length = word_length(w)

    def perm_rep(self, Gamma):
        perms = {g: Gamma.perm(A) for g, A in self.gm.items()}
        for r in self.G.relators():
            for x in range(Gamma.d):
                assert _act_word(perms, r, x) == x, 'relator fails'
        return perms

    def cover(self, Gamma):
        perms = self.perm_rep(Gamma)
        gens = self.G.generators()
        try:
            C = self.M.cover([perms[g] for g in gens])
            conv = 'right'
        except Exception:
            C = self.M.cover([_inv(perms[g]) for g in gens])
            conv = 'left'
        per = self.G.peripheral_curves()
        # orbits: K-cusps <-> cycles of lambda_K (mu_K acts trivially)
        assert all(_act_word(perms, per[self.cuspK][0], x) == x for x in range(Gamma.d))
        lam = [_act_word(perms, per[self.cuspK][1], x) for x in range(Gamma.d)]
        Kcyc = _cycles(lam)
        muT = [_act_word(perms, per[self.cuspT][0], x) for x in range(Gamma.d)]
        Tcyc = _cycles(muT)
        return CoverData(self, Gamma, C, Kcyc, Tcyc, conv)


class CoverData:
    def __init__(self, ml, Gamma, C, Kcyc, Tcyc, conv):
        self.ml, self.Gamma, self.C = ml, Gamma, C
        self.Kcyc, self.Tcyc, self.conv = Kcyc, Tcyc, conv
        self.identify_cusps()

    def identify_cusps(self):
        """Which cusps of C lie over K, their degree m (geodesic = w^m) and the
        meridian slope in SnapPy's (arbitrary, shortest) basis for that cusp.

        SnapPy installs shortest peripheral bases on covers, so the lifted meridian
        is recovered geometrically: a K-cusp of degree m has lattice
        <mu_K, lambda_K^m>, i.e. shape m*s_K up to SL2(Z).  If gamma.(m s_K) = tau_C
        then the meridian is (p,q) = (D,C) where gamma^{-1} = [[A,B],[C,D]]."""
        C, M = self.ml.M, self.C
        sK = complex(self.ml.M.cusp_info(self.ml.cuspK)['shape'])
        shapes = [complex(c['shape']) for c in self.C.cusp_info()]
        nK, nT = len(self.Kcyc), len(self.Tcyc)
        assert len(shapes) == nK + nT, (len(shapes), nK, nT)
        degs = sorted(len(c) for c in self.Kcyc)
        self.shapes = shapes
        # Which cover cusps lie over K, and with what degree: SnapPy builds the cover
        # with tetrahedron i lying over base tetrahedron i // d (checked: this gives a
        # consistent K/T split), so count tetrahedron-vertices at each cusp.
        d = self.Gamma.d
        bidx = self.ml.M._get_cusp_indices_and_peripheral_curve_data()[0]
        cidx = self.C._get_cusp_indices_and_peripheral_curve_data()[0]
        n = len(bidx)
        assert len(cidx) == n * d
        base_cnt = {}
        for t in range(n):
            for v in range(4):
                base_cnt[bidx[t][v]] = base_cnt.get(bidx[t][v], 0) + 1
        over, cnt = {}, {}
        for i in range(len(cidx)):
            for v in range(4):
                c = cidx[i][v]
                bc = bidx[i // d][v]
                if over.setdefault(c, bc) != bc:
                    raise RuntimeError('tetrahedron correspondence inconsistent')
                cnt[c] = cnt.get(c, 0) + 1
        # Meridians: lift mu_K through the covering (exact; see meridian_lift.py).
        # Cusp-shape matching against m*s_K is kept only as an independent cross-check.
        import meridian_lift
        Kdeg = {}
        for i in range(len(shapes)):
            if over[i] == self.ml.cuspK:
                m, r = divmod(cnt[i], base_cnt[self.ml.cuspK])
                assert r == 0
                Kdeg[i] = m
        hint = {}
        for i, m in Kdeg.items():
            cands = _match_slopes(m * sK, shapes[i])
            if cands:
                hint[i] = cands[0]
        lifted, diag = meridian_lift.lifted_meridian_slopes(self.ml.M, self.ml.cuspK, self.C, d, Kdeg, hint)
        Kc = []
        self.slope_checks = {}
        for i, m in sorted(Kdeg.items()):
            sl = lifted[i]
            shape_cands = _match_slopes(m * sK, shapes[i])
            self.slope_checks[i] = dict(lift_integral=diag[i]['ok'],
                                        agrees_with_shape=(sl in shape_cands))
            if sl is None:
                raise RuntimeError(f'cusp {i}: lifted meridian not integral')
            Kc.append((i, m, sl))
        if sorted(m for _, m, _ in Kc) != degs:
            raise RuntimeError(f'cusp identification failed: found {Kc}, expected degrees {degs}')
        self.Kcusps = Kc
        self.meridians_ok = all(v['lift_integral'] and v['agrees_with_shape'] for v in self.slope_checks.values())

    def _resolve_ambiguous(self, Kc, ambiguous):
        """Symmetric cusp lattice: several slopes look like the meridian.  Keep the
        choices for which filling every K-lift gives H1 = Z^(2 + d/6) (= H1 of UT(S))."""
        import itertools
        target = 2 + self.Gamma.d // 6
        good = []
        amb_idx = [a[0] for a in ambiguous]
        for choice in itertools.product(*[a[2] for a in ambiguous]):
            sl = dict(zip(amb_idx, choice))
            N = self.C.copy()
            for (j, m, s0) in Kc:
                N.dehn_fill(sl.get(j, s0), j)
            H = N.homology()
            if H.betti_number() == target and len(H.elementary_divisors()) == target                     and all(e == 0 for e in H.elementary_divisors()):
                good.append(sl)
        if len(good) != 1:
            raise RuntimeError(f'could not resolve meridians: {ambiguous}, {len(good)} candidates pass')
        self.resolved = good[0]
        return [(j, m, good[0].get(j, s0)) for (j, m, s0) in Kc]

    def full_preimage(self):
        C = self.C.copy()
        return C

    def single_lifts(self):
        """For each K-cusp: fill every other K-cusp by its meridian (1,0)."""
        out = []
        for (i, m, _) in self.Kcusps:
            N = self.C.copy()
            for (j, _, sl) in self.Kcusps:
                if j != i:
                    N.dehn_fill(sl, j)
            out.append((i, m, N))
        return out

    def keep(self, keep_idx):
        """Fill every K-cusp not in keep_idx along its meridian."""
        N = self.C.copy()
        for (j, _, sl) in self.Kcusps:
            if j not in keep_idx:
                N.dehn_fill(sl, j)
        return N

    def all_filled(self):
        N = self.C.copy()
        fill = [(0, 0)] * N.num_cusps()      # one dehn_fill call: SnapPea re-solves after every call
        for (j, _, sl) in self.Kcusps:
            fill[j] = sl
        N.dehn_fill(fill)
        return N


def _reduce(tau):
    """Return (tau_red, g) with g in SL2(Z) (as (A,B,C,D)), g.tau = tau_red in std domain."""
    A, B, Cc, D = 1, 0, 0, 1
    for _ in range(10000):
        k = math.floor(tau.real + 0.5)
        if k:
            tau -= k
            A, B = A - k * Cc, B - k * D
        if abs(tau) < 1 - 1e-12:
            tau = -1 / tau
            A, B, Cc, D = -Cc, -D, A, B
        else:
            break
    return tau, (A, B, Cc, D)


def _match_slopes(target, tau, tol=1e-7):
    """All meridian slopes (p,q) (up to sign) such that some gamma in SL2(Z) maps
    target to tau mod Z, i.e. target is the shape in basis (mu'=p mu + q lam, ...)."""
    t1, g1 = _reduce(target)   # g1.target = t1
    t2, g2 = _reduce(tau)      # g2.tau = t2
    out = set()
    # stabiliser/boundary candidates h with h.t1 = t2
    cands = []
    for h in [(1, 0, 0, 1), (1, 1, 0, 1), (1, -1, 0, 1), (0, -1, 1, 0), (1, -1, 1, 0), (0, -1, 1, 1),
              (1, 1, -1, 0), (0, 1, -1, 1), (1, 0, 1, 1), (1, 0, -1, 1)]:
        A, B, Cc, D = h
        if abs(Cc * t1 + D) < 1e-12:
            continue
        v = (A * t1 + B) / (Cc * t1 + D)
        if abs(v - t2) < tol:
            cands.append(h)
    for h in cands:
        # gamma = g2^{-1} h g1 maps target -> tau ; we need gamma^{-1} = g1^{-1} h^{-1} g2
        def mul(P, Q):
            return (P[0] * Q[0] + P[1] * Q[2], P[0] * Q[1] + P[1] * Q[3], P[2] * Q[0] + P[3] * Q[2], P[2] * Q[1] + P[3] * Q[3])
        def inv(P):
            return (P[3], -P[1], -P[2], P[0])
        gi = mul(mul(inv(g1), inv(h)), g2)
        A, B, Cc, D = gi
        p, q = D, Cc
        if p < 0 or (p == 0 and q < 0):
            p, q = -p, -q
        out.add((p, q))
    return sorted(out)
