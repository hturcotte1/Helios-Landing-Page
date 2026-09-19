"""Runs parts (5), (6), (7) of structure_verify.py (parts (1)-(2) ran separately; (3)-(4) superseded by
structure_recover_fast.py / structure_recover_small.py)."""
import structure_verify as V
V.check_pf_expansion(12)
for n in (7, 8, 9, 10): V.check_bilinear(n)
V.check_psi(10)
