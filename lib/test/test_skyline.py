"""SkyLine tests.

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

"""
from __future__ import absolute_import

import pyfits
from numpy.testing import assert_almost_equal, assert_array_less

from .. import skyline

from .test_util import ROOT_DIR
from .test_shared import resolve_imagename

def test_union_simple():
    # Two similar exposures with slight offset and some rotation.
    im1 = resolve_imagename(ROOT_DIR, '2chipA.fits')
    im2 = resolve_imagename(ROOT_DIR, '2chipB.fits')
    
    skyline1 = skyline.SkyLine(im1)
    skyline2 = skyline.SkyLine(im2)
    
    union_1_2 = skyline1.union(skyline2)
    union_2_1 = skyline2.union(skyline1)

    assert_almost_equal(union_1_2.area(), union_2_1.area())

    for m in union_1_2.members:
        assert m in union_2_1.members

    assert len(union_1_2.members) == 4
