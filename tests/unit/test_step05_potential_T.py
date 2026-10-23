"""step05 potential temperature (UNESCO 1983) against published check values."""
import numpy as np

import step05


def test_adiabatic_gradient_unesco_check_value():
    # UNESCO 1983 (Fofonoff & Millard) check value: S=40, T=40 C, P=10000 dbar -> 3.255976e-4 C/dbar
    assert abs(step05.calculate_adiabatic_T_gradient(40.0, 40.0, 10000.0) - 3.255976e-4) < 1e-9
