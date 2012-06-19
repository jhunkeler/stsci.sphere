"""
SkyLine tests.

:Author: Pey Lian Lim

:Organization: Space Telescope Science Institute

Examples
--------
>>> cd path/to/project
>>> nosetests

"""
from __future__ import absolute_import

from copy import copy
from numpy.testing import assert_almost_equal

from ..skyline import SkyLine

from .test_util import ROOT_DIR
from .test_shared import resolve_imagename


#---------------------------------------#
# Load footprints used by all the tests #
#---------------------------------------#
f_2chipA  = resolve_imagename(ROOT_DIR, '2chipA.fits') # ACS/WFC #1
im_2chipA = SkyLine(f_2chipA)
f_2chipB  = resolve_imagename(ROOT_DIR, '2chipB.fits') # ACS/WFC #2
im_2chipB = SkyLine(f_2chipB)
f_2chipC  = resolve_imagename(ROOT_DIR, '2chipC.fits') # WFC3/UVIS
im_2chipC = SkyLine(f_2chipC)
f_66_tan  = resolve_imagename(ROOT_DIR, '1904-66_TAN.fits')
im_66_tan = SkyLine(f_66_tan, extname='primary')


#----- SHARED FUNCTIONS -----

def compare_members(im1, im2):
    assert len(im1.members) == len(im2.members)

    for m in im1.members:
        assert m in im2.members

    for m in im2.members:
        assert m in im1.members


#----- MEMBERSHIP -----

def do_member_overlap(im):
    for m in im.members:
        assert_almost_equal(m.polygon.overlap(im), 1.0)

def test_membership():
    do_member_overlap(im_2chipA)
    do_member_overlap(im_2chipB)
    do_member_overlap(im_66_tan)

    assert len(im_2chipA.members) == 2
    assert im_2chipA.members[0].fname == f_2chipA
    assert im_2chipA.members[0].ext == 1
    assert im_2chipA.members[1].fname == f_2chipA
    assert im_2chipA.members[1].ext == 4


#----- COPY -----

def test_copy():
    a_copy = copy(im_2chipA)
    assert a_copy is not im_2chipA


#----- SPHERICAL POLYGON RELATED -----

def test_sphericalpolygon():
    assert im_2chipA.contains_point(im_2chipA.inside)
    assert im_2chipA.intersects_poly(im_2chipB.polygon)
    assert im_2chipA.intersects_arc(im_2chipA.inside, im_2chipB.inside)
    assert im_2chipA.overlap(im_2chipB) < im_2chipA.overlap(im_2chipA)
    assert_almost_equal(im_2chipA.area(), im_2chipB.area())


#----- UNION -----

def do_add_image(im1, im2):
    u1 = im1.add_image(im2)
    u2 = im2.add_image(im1)

    assert u1.same_points_as(u2)
    compare_members(u1, u2)

def test_add_image():    
    # Dithered
    do_add_image(im_2chipA, im_2chipB)

    # Not related
    do_add_image(im_2chipA, im_66_tan)


# ----- INTERSECTION -----

def do_intersect_image(im1, im2):
    i1 = im1.find_intersection(im2)
    i2 = im2.find_intersection(im1)

    assert i1.same_points_as(i2)
    compare_members(i1, i2)

def test_find_intersection():   
    # Dithered
    do_intersect_image(im_2chipA, im_2chipB)

    # Not related
    do_intersect_image(im_2chipA, im_66_tan)


# ----- INTENDED USE CASE -----

def test_science():
    skylines = [im_2chipA, im_2chipB, im_2chipC, im_66_tan]

    # TODO: Add Warren's example use case


# ----- UNSTABLE -----

def DISABLED_unstable_overlap():
    i1 = im_2chipA.find_intersection(im_2chipB)
    i2 = im_2chipB.find_intersection(im_2chipA)
    
    u1 = im_2chipA.add_image(im_2chipB)
    u2 = im_2chipB.add_image(im_2chipA)
    
    assert_almost_equal(i1.overlap(u1), 1.0)
    assert_almost_equal(i1.overlap(i2), 1.0) # failed here - known bug
    assert_almost_equal(u1.overlap(u2), 1.0)
