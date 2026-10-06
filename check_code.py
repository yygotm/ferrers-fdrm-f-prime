"""Independent checker for an [F', 6, 7]_q FDRM code, F' = [1,2,3,4,8,8,8,8,8].

Input file: 6 non-empty lines, one basis matrix per line, 8 rows separated by spaces,
each row 9 symbols in 0..q-1 (column 0 first). For q = 4 the symbols are
0, 1, 2 = w, 3 = w^2 = w + 1 in GF(4) = GF(2)[w]/(w^2 + w + 1).

Checks: (1) every basis matrix is supported on F' (top-aligned columns of heights C);
(2) the 6 basis matrices are linearly independent; (3) every nonzero F_q-combination
(all q^6 - 1 of them, no projective reduction) has rank >= 7.
usage: python check_code.py q file
"""
import sys
from itertools import product

C = [1, 2, 3, 4, 8, 8, 8, 8, 8]
NR, NC, K, D = 8, 9, 6, 7


def field(q):
    if q in (2, 3, 5, 7):
        add = [[(a + b) % q for b in range(q)] for a in range(q)]
        mul = [[(a * b) % q for b in range(q)] for a in range(q)]
    elif q == 4:
        # elements as bit pairs (c1 c0) meaning c1*w + c0; 2 = w, 3 = w + 1 = w^2
        def m(a, b):
            r = 0
            for i in range(2):
                if b >> i & 1:
                    r ^= a << i
            if r & 4:
                r ^= 0b111  # w^2 = w + 1
            return r
        add = [[a ^ b for b in range(4)] for a in range(4)]
        mul = [[m(a, b) for b in range(4)] for a in range(4)]
    else:
        raise SystemExit("unsupported q")
    inv = [None] + [next(b for b in range(1, q) if mul[a][b] == 1) for a in range(1, q)]
    neg = [next(b for b in range(q) if add[a][b] == 0) for a in range(q)]
    return add, mul, inv, neg


def rank(rows, q, F):
    add, mul, inv, neg = F
    a = [list(r) for r in rows]
    r = 0
    ncol = len(a[0])
    for c in range(ncol):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        s = inv[a[r][c]]
        a[r] = [mul[s][x] for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                f = neg[a[i][c]]
                a[i] = [add[x][mul[f][y]] for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def main():
    q = int(sys.argv[1])
    F = field(q)
    add, mul, _, _ = F
    lines = [l.split() for l in open(sys.argv[2]) if l.strip()]
    assert len(lines) == K, f"expected {K} basis matrices, got {len(lines)}"
    B = []
    for l in lines:
        assert len(l) == NR and all(len(r) == NC for r in l)
        M = [[int(ch) for ch in r] for r in l]
        assert all(0 <= x < q for row in M for x in row)
        for i in range(NR):
            for j in range(NC):
                if M[i][j] and i >= C[j]:
                    raise SystemExit(f"FAIL: entry ({i},{j}) outside F'")
        B.append(M)
    flat = [[x for row in M for x in row] for M in B]
    rk = rank(flat, q, F)
    print("basis rank", rk)
    assert rk == K, "FAIL: basis dependent"
    hist = {}
    for coef in product(range(q), repeat=K):
        if not any(coef):
            continue
        M = [[0] * NC for _ in range(NR)]
        for a, Bk in zip(coef, B):
            if a:
                for i in range(NR):
                    for j in range(NC):
                        M[i][j] = add[M[i][j]][mul[a][Bk[i][j]]]
        r = rank(M, q, F)
        hist[r] = hist.get(r, 0) + 1
    print("rank distribution over", q ** K - 1, "nonzero codewords:", dict(sorted(hist.items())))
    ok = min(hist) >= D
    print("PASS" if ok else "FAIL", f"[F',{K},{D}]_{q}")
    sys.exit(0 if ok else 1)


main()
