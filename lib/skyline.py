# -*- coding: utf-8 -*-

# Copyright (C) 2011 Association of Universities for Research in
# Astronomy (AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#     1. Redistributions of source code must retain the above
#       copyright notice, this list of conditions and the following
#       disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Manage outlines on the sky.

This module provides support for working with footprints
on the sky. Primary use case would use the following
generalized steps::

    #. Initialize `SkyLine` objects for each input image.
       This object would be the union of all the input
       image's individual chips WCS footprints.

    #. Determine overlap between all images. The
       determination would employ a recursive operation
       to return the extended list of all overlap values
       computed as [img1 vs [img2,img3,...,imgN],img2 vs
       [img3,...,imgN],...]

    #. Select the pair with the largest overlap, or the
       pair which produces the largest overlap with the
       first input image. This defines the initial
       reference `SkyLine` object.

    #. Perform some operation on the 2 images: for example,
       match sky in intersecting regions, or aligning
       second image with the first (reference) image.

    #. Update the second image, either apply the sky value
       or correct the WCS, then generate a new `SkyLine`
       object for that image.

    #. Create a new reference `SkyLine` object as the union
       of the initial reference object and the newly
       updated `SkyLine` object.

    #. Repeat Steps 2-6 for all remaining input images.

This process will work reasonably fast as most operations
are performed using the `SkyLine` objects and WCS information
solely, not image data itself.

.. note:: Requires Python 2.7 or later.

:Authors: Pey Lian Lim, Warren Hack, Michael Droettboom

:Organization: Space Telescope Science Institute

:History:
    * 2012-05-25 PLL started coding. Original class
      structure by WJH. Parent class by MD.

"""
from __future__ import division, print_function, absolute_import

# STDLIB
from copy import copy, deepcopy

# THIRD-PARTY
import pyfits
from stwcs import wcsutil

# LOCAL
from .polygon import SphericalPolygon

__all__ = ['SkyLine']
__version__ = '0.2a'
__vdate__ = '12-Jun-2012'

class SkyLineMember(object):

    def __init__(self, fname, ext):
        """
        Container for `SkyLine` members.

        Given FITS image and extension, will store its
        `SphericalPolygon` instance from WCS under
        `polygon`.

        Parameters
        ----------
        self: obj
            SkyLineMember instance.

        fname: str
            FITS image.

        ext: int
            Image extension.

        """
        self._fname = fname
        self._ext = ext
        self._polygon = SphericalPolygon.from_wcs(
            wcsutil.HSTWCS(fname, ext=ext))

    def __repr__(self):
        return 'SkyLineMember(%r, %r, %r)' % (self.fname, self.ext,
                                              self.polygon)

    @property
    def fname(self):
        return self._fname

    @property
    def ext(self):
        return self._ext

    @property
    def polygon(self):
        return self._polygon

class SkyLine(object):

    def __init__(self, fname, extname='SCI'):
        """
        Initialize `SkyLine` object instance.

        Parameters
        ----------
        self: obj
            `SkyLine` instance.

        fname: str
            FITS image. `None` to create empty `SkyLine`.

        extname: str
            EXTNAME to use. SCI is recommended for normal
            HST images. PRIMARY if image is single ext.

        """
        extname = extname.upper()
        
        # Convert SCI data to SkyLineMember
        if fname is not None:
            with pyfits.open(fname) as pf:
                self.members = [SkyLineMember(fname, i)
                                for i,ext in enumerate(pf)
                                if extname in ext.name.upper()]
        else:
            self.members = []

        # Put mosaic of all the chips in SkyLine
        if len(self.members) > 0:
            self.polygon = SphericalPolygon.multi_union([m.polygon for m in self.members])
        else:
            self.polygon = SphericalPolygon([])

    def __getattr__(self, what):
        """Control attribute access to `SphericalPolygon`."""
        if what in ('from_radec', 'from_cone', 'from_wcs',
                    'multi_union', 'multi_intersection',
                    '_find_new_inside',):
            raise AttributeError('\'SkyLine\' object has no attribute \'%s\'' %
                                 what)
        else:
            return getattr(self.polygon, what)

    def __copy__(self):
        return deepcopy(self)
    
    def __repr__(self):
        return 'SkyLine(%r, %r)' % (self.polygon, self.members)

    @property
    def polygon(self):
        """`SphericalPolygon` portion of `SkyLine`."""
        return self._polygon

    @polygon.setter
    def polygon(self, value):
        """Deep copy a `SphericalPolygon`."""
        assert isinstance(value, SphericalPolygon)
        self._polygon = copy(value)

    @property
    def members(self):
        """List of `SkyLineMember` objects."""
        return self._members

    @members.setter
    def members(self, values):
        """Make sure `SkyLineMember` entries are unique."""
        self._members = []

        # Not using set to preserve order
        for v in values:
            # Report corrupted members list instead of skipping
            assert isinstance(v, SkyLineMember)

            if v not in self._members:
                self._members.append(v)

    def _find_members(self, given_members):
        """
        Find `SkyLineMember` in *given_members* that is in
        *self*. This is used for intersection.

        Parameters
        ----------
        self: obj
            `SkyLine` instance.

        given_members: list
            List of `SkyLineMember` to consider.

        Returns
        -------
        new_members: list
            List of `SkyLineMember` belonging to *self*.

        """
        if len(self.points) > 0:
            out_mem = [m for m in given_members if
                       self.intersects_poly(m.polygon)]
        else:
            out_mem = []
        return out_mem

    def add_image(self, other):
        """
        Return a new `SkyLine` that is the union of *self*
        and *other*.

        .. warning:: `SkyLine.union` only returns `polygon`
            without `members`.

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.add_image(s2)

        """
        newcls = self.__class__(None)
        newcls.polygon = self.union(other)
        newcls.members = self.members + other.members
        return newcls

    def find_intersection(self, other):
        """
        Return a new `SkyLine` that is the intersection of
        *self* and *other*.

        .. warning:: `SkyLine.intersection` only returns
            `polygon` without `members`.

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.find_intersection(s2)

        """
        newcls = self.__class__(None)
        newcls.polygon = self.intersection(other)
        newcls.members = newcls._find_members(self.members + other.members)
        return newcls
