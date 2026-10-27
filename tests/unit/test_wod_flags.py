"""WOD quality flags (ecco_insitu.wod_flags): cast numbers from prof_descr, and flags on our depths."""
import numpy as np

from ecco_insitu import wod_flags as W


def test_cast_numbers():
    codes, casts = W.cast_numbers([b'WOD-XBT//3975702:BS142898_7716', b'WOD-MRB//22870282_52004_Tropic',
                                   'WOD-OSD//13512574:_ ', b'ARGO-PFL-aoml//D5902328-3', b'MEOP-D00050144-6/hd2'])
    assert list(codes) == ['XBT', 'MRB', 'OSD', '', '']
    assert list(casts) == [3975702, 22870282, 13512574, -1, -1]


def test_level_flags():
    z_std = np.array([0., 10., 20., 30., 50.])
    f_std = np.array([0, 0, 3, -1, 0], np.int8)        # 30 m has no value: skipped
    z = np.array([0., 5., 10., 15., 20., 25., 40., 50., 60.])
    #   exact 0 -> 0; 5 between 0,10 -> 0; exact 10 -> 0; 15 between 10,20 -> 3; exact 20 -> 3;
    #   25 and 40 between 20 and 50 (30 skipped) -> 3; exact 50 -> 0; 60 beyond the cast -> 0
    assert list(W.level_flags(z_std, f_std, z)) == [0, 0, 0, 3, 3, 3, 3, 0, 0]


def test_level_flags_no_values():
    assert list(W.level_flags(np.array([0., 10.]), np.array([-1, -1], np.int8), np.array([5.]))) == [0]
