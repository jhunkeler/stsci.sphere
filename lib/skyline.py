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
__version__ = '0.1a'
__vdate__ = '07-Jun-2012'

class SkyLineMember(object):

    def __init__(self, fname, ext):
        """
        Container for `SkyLine` members.

        Given FITS image and extension, will store its
        `SphericalPolygon` instance from WCS under
        `polygon`.

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

class SkyLine(SphericalPolygon):

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
        SphericalPolygon.__init__(self, [])

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
            self.polygon = SphericalPolygon.multi_union(
                [m.polygon for m in self.members])

    def __repr__(self):
        return 'SkyLine(%r, %r, %r)' % (self.points, self.inside, self.members)

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
                       self.contains_point(m.polygon.inside)]
        else:
            out_mem = []
        return out_mem

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

    @property
    def polygon(self):
        """`SphericalPolygon` portion of `SkyLine`."""
        return SphericalPolygon(self.points, self.inside)

    @polygon.setter
    def polygon(self, sph):
        """Set `SkyLine` to given `SphericalPolygon` properties."""
        assert isinstance(sph, SphericalPolygon)
        self._points = sph.points
        self._inside = sph.inside

    @staticmethod
    def _sep_poly_mem(skylines):
        """
        Separate polygons and members of given *skylines*
        for further processing.

        Returns
        -------
        all_poly: list
            List of `SphericalPolygon`.

        all_mem: list
            List of `SkyLineMember`.

        """
        all_poly, all_mem = [], []
        
        for a in skylines:
            all_poly.append(a.polygon)
            all_mem += a.members

        return all_poly, all_mem

    @classmethod
    def _overload_parentcls(cls, mem, func, *args, **kwargs):
        """Call `SphericalPolygon` class method but return `SkyLine`."""
        newcls = cls(None)
        newcls.polygon = func(*args, **kwargs)
        newcls.members = mem
        return newcls

    @classmethod
    def from_radec(cls, *args, **kwargs):
        """
        Create a new `SkyLine` from a list of (*ra*, *dec*)
        points.

        See also
        --------
        sphere.polygon.SphericalPolygon.from_radec

        """
        return cls._overload_parentcls([], SphericalPolygon.from_radec,
                                       *args, **kwargs)

    @classmethod
    def from_cone(cls, *args, **kwargs):
        """
        Create a new `SkyLine` from a cone (otherwise known
        as a 'small circle') defined using (*ra*, *dec*, *radius*).

        See also
        --------
        sphere.polygon.SphericalPolygon.from_cone

        """
        return cls._overload_parentcls([], SphericalPolygon.from_cone,
                                       *args, **kwargs)

    @classmethod
    def from_wcs(cls, *args, **kwargs):
        """
        Create a new `SkyLine` from the footprint of a FITS
        WCS specification.

        See also
        --------
        sphere.polygon.SphericalPolygon.from_wcs

        """
        return cls._overload_parentcls([], SphericalPolygon.from_wcs,
                                       *args, **kwargs)

    @classmethod
    def multi_union(cls, skylines, **kwargs):
        """
        Return a new `SkyLine` that is the union of all of
        the *skylines*.

        See also
        --------
        sphere.polygon.SphericalPolygon.multi_union

        """
        all_poly, all_mem = cls._sep_poly_mem(skylines)
        return cls._overload_parentcls(all_mem, SphericalPolygon.multi_union,
                                       all_poly, **kwargs)

    @classmethod
    def multi_intersection(cls, skylines, **kwargs):
        """
        Return a new `SkyLine` that is the intersection of
        all of the *skylines*.

        See also
        --------
        sphere.polygon.SphericalPolygon.multi_intersection

        """
        all_poly, all_mem = cls._sep_poly_mem(skylines)
        newcls = cls._overload_parentcls([],
                                         SphericalPolygon.multi_intersection,
                                         all_poly, **kwargs)
        newcls.members = newcls._find_members(all_mem)
        return newcls

    def union(self, other):
        """
        Return a new `SkyLine` that is the union of *self*
        and *other*.

        Parameters
        ----------
        self, other: obj
            `SkyLine` instance.

        Returns
        -------
        out_skyline: obj
            `SkyLine` instance.

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.union(s2)

        See also
        --------
        sphere.polygon.SphericalPolygon.union

        """
        out_skyline = self.__class__(None)
        out_skyline.polygon = self.polygon.union(other.polygon)
        out_skyline.members = self.members + other.members
        return out_skyline

    def intersection(self, other):
        """
        Return a new `SkyLine` that is the intersection of
        *self* and *other*.

        Parameters
        ----------
        self, other: obj
            `SkyLine` instance.

        Returns
        -------
        out_skyline: obj
            `SkyLine` instance.

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s3 = s1.intersection(s2)

        See also
        --------
        sphere.polygon.SphericalPolygon.intersection

        """
        out_skyline = self.__class__(None)
        out_skyline.polygon = self.polygon.intersection(other.polygon)
        out_skyline.members = out_skyline._find_members(
            self.members + other.members)
        return out_skyline
