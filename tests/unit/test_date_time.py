"""NCEI._normalize_date_time: out-of-range WOD dates/times roll forward, as prof_date does."""
import numpy as np
import pytest

NCEI = pytest.importorskip('NCEI')


@pytest.mark.parametrize('ymd, hms, want', [
    (20100101, 42760, (20100101, 42800)),      # SS == 60
    (20100101, 46000, (20100101, 50000)),      # MM == 60
    (20100101, 235960, (20100102, 0)),         # 23:59:60 -> next day
    (20100101, 240000, (20100102, 0)),         # 24:00:00 midnight convention
    (20100931, 120000, (20101001, 120000)),    # Sep 31 -> Oct 1
    (20101231, 240000, (20110101, 0)),         # year rollover
    (20120228, 240000, (20120229, 0)),         # leap year
    (20100615, 123456, (20100615, 123456)),    # valid: unchanged
])
def test_rollover(ymd, hms, want):
    d, t = NCEI._normalize_date_time(np.array([ymd]), np.array([hms]))
    assert (int(d[0]), int(t[0])) == want


@pytest.mark.parametrize('ymd, hms', [(-9999, -9999), (20101301, 0), (20100100, 0), (20100101, 250000)])
def test_implausible_values_untouched(ymd, hms):
    d, t = NCEI._normalize_date_time(np.array([ymd]), np.array([hms]))
    assert (int(d[0]), int(t[0])) == (ymd, hms)
