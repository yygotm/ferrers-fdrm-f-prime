# Optimal [F', 6, 7]_q FDRM codes for q = 2, 3, 4, 5

This repository gives explicit optimal Ferrers diagram rank-metric (FDRM) codes for the
Ferrers diagram

    F' = [1, 2, 3, 4, 8, 8, 8, 8, 8]

over F_2, F_3, F_4 and F_5. Together with the construction for q >= 7 in
Fu, Gao, Mesnager and Wang, *Optimal Ferrers Diagram Rank-Metric Codes: New Constructions,
Diagram Combinations, and Applications to Constant-Dimension Subspace Codes*
(arXiv:2609.24239v1, Remark 2), this answers Open Problem 1 of that paper affirmatively
for every prime power q, as far as this computation goes.

## The question

F' has 8 rows and 9 columns. Its column heights, from left to right, are 1, 2, 3, 4, 8, 8, 8,
8, 8, with the dots of each column top-aligned (50 dots). It is obtained from
F = [1, 2, 3, 4, 8, 8, 8, 8, 9] by deleting one "vain" dot. An [F', 6, 7]_q code is a
6-dimensional F_q-linear space of 8 x 9 matrices supported on the dots of F', in which every
nonzero matrix has rank at least 7.

The Etzion–Silberstein bound gives

    v_min(F', 7) = min_{0<=i<=6} v_i = min(6, 6, 9, 11, 12, 12, 10) = 6,

where v_i is the number of dots of F' left after deleting the top i rows and the rightmost
6 - i columns. So a code of dimension 6 is optimal.

## Result

| q | basis file | nonzero codewords | rank 7 | rank 8 |
|--:|:--|--:|--:|--:|
| 2 | `codes/code_q2.txt` | 63 | 50 | 13 |
| 3 | `codes/code_q3.txt` | 728 | 350 | 378 |
| 4 | `codes/code_q4.txt` | 4095 | 1140 | 2955 |
| 5 | `codes/code_q5.txt` | 15624 | 3048 | 12576 |

Each file has 6 lines, one basis matrix per line. A line consists of 8 words (the rows), and
the j-th character of a word is the entry in column j (column 0 first). For q = 2, 3, 5 the
characters 0..q-1 are the elements of F_q = Z/qZ. For q = 4 the characters 0, 1, 2, 3 mean
0, 1, w, w^2 in F_4 = F_2[w]/(w^2 + w + 1), so w^2 = w + 1. SHA-256 sums are in
`codes/SHA256SUMS`.

## Checking

```
python check_code.py 2 codes/code_q2.txt
python check_code.py 3 codes/code_q3.txt
python check_code.py 4 codes/code_q4.txt
python check_code.py 5 codes/code_q5.txt
```

or `sh run_all.sh`. `check_code.py` (Python 3, standard library only) checks that

1. every basis matrix is supported on F',
2. the 6 basis matrices are linearly independent over F_q,
3. all q^6 - 1 nonzero F_q-combinations (no projective reduction) have rank >= 7,

and prints the rank distribution. The outputs are in `results/`. The whole check takes a few
seconds.

The four codes were also re-checked by a separately written program (another AI session,
which was not shown `check_code.py` or the search code). It verified the same three
properties twice, by Gaussian elimination and by explicit minors (a nonzero 7 x 7 minor for
every nonzero codeword), found the same rank distributions, and recomputed the bound
v_min(F', 7) = 6. Sage or Magma was not used.

## How the codes were found

By the bound at i = 0, the projection of the code onto the 6 dots of columns 0..2 must be
injective. So the code can be written as the graph of a linear map from those 6 dots to the
other 44 dots: basis matrix k has a 1 at the k-th dot of columns 0..2, zeros at the other
five of those dots, and free entries elsewhere. The free entries were found by simulated
annealing that minimises the sum over nonzero codewords (up to scalars) of max(0, 7 - rank).

- `search/saq.c`: general q in {2, 3, 4, 5, 7}. Usage: `saq q seed steps T0`.
- `search/sa2.c`: the earlier binary version. Usage: `sa2 seed steps T0`.

```
gcc -O2 -o saq search/saq.c -lm
./saq 5 1 2000000 2
```

The parameters used for the published codes were not recorded, so a rerun gives different
codes. With `seed = 1, steps = 2000000, T0 = 2` (gcc 15.2, Ubuntu), the search finds a code
for each q = 2, 3, 4, 5 within 4 to 15 seconds, and each passes `check_code.py`.

## Scope

This is a computer search with an exhaustive check of four explicit codes over four fields.
It does not give a uniform construction for all q. The computations were carried out with
the help of an AI assistant (Anthropic's Claude) and checked independently with OpenAI's
GPT. It is not a formally verified proof.

## Author

Hirotaka Shimizu
