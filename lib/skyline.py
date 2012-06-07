"""Manage outlines on the sky.

This module provides support for working with footprints
on the sky. Primary use case would use the following
generalized steps::

    #. Initialize SkyLine objects for each input image.
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
       reference SkyLine object.

    #. Perform some operation on the 2 images: for example,
       match sky in intersecting regions, or aligning
       second image with the first (reference) image.

    #. Update the second image, either apply the sky value
       or correct the WCS, then generate a new SkyLine
       object for that image.

    #. Create a new reference SkyLine object as the union
       of the initial reference object and the newly
       updated SkyLine object.

    #. Repeat Steps 2-6 for all remaining input images.

This process will work reasonably fast as most operations
are performed using the SkyLine objects and WCS information
solely, not image data itself.

.. note:: Requires Python 2.7 or later.

:Authors: Pey Lian Lim, W. Hack, M. Droettboom

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
__vdate__ = '06-Jun-2012'

class SkyLineMember(object):

    def __init__(self, fname, ext):
        """Container for SkyLine members.

        Given FITS image and extension, will store its SphericalPolygon
        instance from WCS under `polygon`.

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
        """Initialize SkyLine object instance.

        Parameters
        ----------
        self: obj
            SkyLine instance.

        fname: str
            FITS image. None to create empty SkyLine.

        extname: str
            EXTNAME to use. SCI is recommended for normal
            HST images. PRIMARY if image is single ext.

        """
        SphericalPolygon.__init__(self, [])

        # Convert SCI data to SkyLineMember
        if fname is not None:
            with pyfits.open(fname) as pf:
                new_members = [SkyLineMember(fname, i)
                               for i,ext in enumerate(pf)
                               if extname in ext.name.upper()]
        else:
            new_members = []

        # Put mosaic of all the chips in SkyLine
        if len(new_members) > 0:
            new_polygon = SphericalPolygon.multi_union(
                [m.polygon for m in new_members])
        # Empty class
        else:
            new_polygon = self

        self._update(new_polygon, new_members)

    def __repr__(self):
        return 'SkyLine(%r, %r, %r)' % (self.points, self.inside, self.members)

    def _update(self, new_polygon, new_members):
        """
        Update *self* attributes to use given polygon and
        new members.

        Parameters
        ----------
        self: obj
            SkyLine instance to update.

        new_polygon: obj
            SphericalPolygon instance to use.

        new_members: list
            List of SkyLineMember associated with `new_polygon`.

        """
        self._points = new_polygon.points
        self._inside = new_polygon.inside
        self._members = new_members

    def _find_new_members(self, other):
        """
        Find SkyLineMember that is in *other* but not in *self*.

        This is used internally to make sure there are no duplicate
        SkyLineMember entries. Order is preserved, with *self*
        listed first, followed by each new member from *other*.

        Parameters
        ----------
        self, other: obj
            `SkyLine` instance.

        Returns
        -------
        List of SkyLineMember that qualifies.

        """
        return [m for m in other.members if m not in self.members]

    @property
    def members(self):
        """List of SkyLineMember objects."""
        return self._members

    @property
    def polygons(self):
        """List of SkyLineMember polygons."""
        return [m.polygon for m in self.members]

    @property
    def polygon(self):
        """SphericalPolygon portion of SkyLine."""
        return SphericalPolygon(self.points, self.inside)

    @classmethod
    def _overload_parentcls(cls, func, *args, **kwargs):
        """Call SphericalPolygon class method but return SkyLine."""
        newcls = cls(None)
        newcls._update(func(*args, **kwargs), None)
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
        return cls._overload_parentcls(SphericalPolygon.from_radec,
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
        return cls._overload_parentcls(SphericalPolygon.from_cone,
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
        return cls._overload_parentcls(SphericalPolygon.from_wcs,
                                       *args, **kwargs)

    @classmethod
    def multi_union(cls, *args, **kwargs):
        """
        Return a new `SkyLine` that is the union of all of the
        polygons in *polygons*.

        See also
        --------
        sphere.polygon.SphericalPolygon.multi_union

        """
        return cls._overload_parentcls(SphericalPolygon.multi_union,
                                       *args, **kwargs)

    @classmethod
    def multi_intersection(cls, *args, **kwargs):
        """
        Return a new `SkyLine` that is the intersection of
        all of the polygons in *polygons*.

        See also
        --------
        sphere.polygon.SphericalPolygon.multi_intersection

        """
        return cls._overload_parentcls(SphericalPolygon.multi_intersection,
                                       *args, **kwargs)

    def union(self, other):
        """
        Return a new `SkyLine` that is the union of *self* and *other*.
        Skips *other* members that are already in *self*.

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
        out_skyline = copy(self)
        new_members = self._find_new_members(other)

        if len(new_members) > 0:
            out_skyline._update(self.polygon.union(other.polygon),
                                self.members + new_members)

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
        new_members = self._find_new_members(other)
        out_sph = self.polygon.intersection(other.polygon)

        if len(out_sph.points) > 0:
            new_members = [m for m in (self.members + new_members) if
                           out_sph.contains_point(m.polygon.inside)]
            out_skyline._update(out_sph, new_members)

        return out_skyline
