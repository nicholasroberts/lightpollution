"""Numerical audit of recovered R response equations, not site-PRC reproduction.

Standard-library only. Synthetic inputs; no field images are processed.
Run: python src/audit_original_spider_response.py --self-test
"""

import argparse
import hashlib
import json
import math
from pathlib import Path


def literal_response(p0, p1):
    """Literal scalar translation of the supplied R par/perp expressions."""
    r = 10 ** (-0.4)
    q = 10 ** (-0.04)
    perp = p0 - 0.5 * (q * (p0 - p1) + r * (p0 + p1))
    par = p0 - 0.5 * (q * (p0 + p1) + r * (p0 - p1))
    return par, perp


def self_test():
    q, r = 10 ** (-0.04), 10 ** (-0.4)
    max_error = 0.0
    for p0 in (0.0, 0.1, 1.0, 10.0):
        for fraction in (-1.0, -0.5, 0.0, 0.5, 1.0):
            p1 = p0 * fraction
            par, perp = literal_response(p0, p1)
            expected_diff = -(q - r) * p1
            expected_sum = (2 - q - r) * p0
            max_error = max(max_error, abs(par - perp - expected_diff),
                            abs(par + perp - expected_sum))
            assert math.isclose(par - perp, expected_diff, abs_tol=1e-13)
            assert math.isclose(par + perp, expected_sum, abs_tol=1e-13)
            assert par >= -1e-13 and perp >= -1e-13
            # Unpolarized added intensity changes P0 but not P1.
            added_par, added_perp = literal_response(p0 + 2.0, p1)
            assert math.isclose(added_par - added_perp, par - perp,
                                abs_tol=1e-13)
    par, perp = literal_response(1.0, 0.5)
    swapped_par, swapped_perp = literal_response(1.0, -0.5)
    assert math.isclose(par, swapped_perp, abs_tol=1e-13)
    assert math.isclose(perp, swapped_par, abs_tol=1e-13)
    source = (Path(__file__).resolve().parents[1] / 'archive/legacy_prc/code/'
              'sky_image_processing_supplied_2026_11_05.R')
    return {
        'status': 'PASS',
        'scope': 'synthetic algebra audit; R execution and field reproduction NOT performed',
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'a_par': 0.4, 'a_perp': 0.04, 'sent_r': r, 'sent_q': q,
        'difference_coefficient': -(q-r),
        'sum_coefficient': 2-q-r,
        'maximum_absolute_identity_error': max_error,
        'example_P0_1_P1_0p5': {'par': par, 'perp': perp,
                                  'pol_im': par-perp, 'int_im': par+perp},
        'unpolarized_addition': 'raw pol.im unchanged; final eye-to-eye metric unrecovered',
        'site_PRC_reproduced': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(self_test(), indent=2))
