"""step11 spatial-scaling time windows."""
import pytest

import step11_spatial_scaling as s11


@pytest.mark.parametrize('ymd, want', [(20150101, (2015, 1)), (20151231, (2015, 12)), (19920615, (1992, 6))])
def test_monthly(ymd, want):
    assert s11._window_monthly(ymd) == want


@pytest.mark.parametrize('ymd, want', [
    (20150101, (2015, 0)), (20150110, (2015, 0)), (20150111, (2015, 1)),
    (20151231, (2015, 36)), (20161231, (2016, 36)),
])
def test_10day(ymd, want):
    assert s11._window_10day(ymd) == want
