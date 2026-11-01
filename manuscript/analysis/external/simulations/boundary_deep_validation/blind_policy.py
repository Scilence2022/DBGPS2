"""Locked policies. This module never receives a reference, oracle, or S/W label."""
import numpy as np


def select(status, candidates, settings):
    assert status.ndim == 2 and candidates.shape[:2] == status.shape
    ok = status == 0
    counts = ok.sum(axis=1)
    rkey = lambda j: settings[j]['r'] if settings[j]['r'] > 0 else float('inf')
    c1 = [j for j, s in enumerate(settings) if s['c'] == 1]
    default = next(j for j in c1 if settings[j]['r'] == 2.55)
    best_r = min(c1, key=lambda j: (-int(counts[j]), rkey(j)))
    best_cr = min(range(len(settings)), key=lambda j: (-int(counts[j]), -settings[j]['c'], rkey(j)))
    out = {}
    for name, j in [('default', default), ('global_crc_r', best_r), ('global_crc_cr', best_cr)]:
        out[name] = dict(accepted=ok[j].copy(), candidates=candidates[j].copy(), selected_combo=j)
    chosen = np.zeros_like(candidates[0]); seen = np.zeros(status.shape[1], dtype=np.int16)
    conflict = np.zeros(status.shape[1], dtype=bool)
    for j in c1:
        accepted = ok[j]
        conflict |= accepted & (seen > 0) & np.any(candidates[j] != chosen, axis=1)
        first = accepted & (seen == 0)
        chosen[first] = candidates[j, first]
        seen += accepted
    for name, support in [('union_consensus', 1), ('two_support_consensus', 2)]:
        out[name] = dict(accepted=(seen >= support) & ~conflict,
                         candidates=chosen.copy(), selected_combo=None,
                         conflicts=int(conflict.sum()))
    return out
