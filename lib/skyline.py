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

:Authors: Pey Lian Lim, W. Hack

:Organization: Space Telescope Science Institute

:History:
    * 2012-05-25 PLL updated doc. Original class structure by WJH.

Examples
--------
>>> from sphere import SkyLine

"""
from __future__ import division, print_function

# STDLIB
from copy import deepcopy

# THIRD-PARTY
import pyfits
from stwcs import wcsutil

# LOCAL
from sphere.polygon import SphericalPolygon

__all__ = ['SkyLine']
__version__ = '0.1a'
__vdate__ = '31-May-2012'

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

    def __init__(self, fname):
        """Initialize SkyLine object instance.

        Parameters
        ----------
        self: obj
            SkyLine instance.

        fname: str
            FITS image.

        """
        # Inherit from SphericalPolygon
        SphericalPolygon.__init__(self, [], None)
        
        # Convert SCI data to SkyLineMember
        poly_list = []
        with pyfits.open(fname) as pf:
            for i,ext in enumerate(pf):
                if 'SCI' in ext.name.upper():
                     poly_list.append(SkyLineMember(fname, i))

        assert len(poly_list) > 0, \
            'SkyLine cannot find SCI ext in {}.'.format(fname)

        # Put mosaic of all the chips in SkyLine
        self._update(self.multi_union([m.polygon for m in poly_list]),
                     poly_list)

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

    @property
    def members(self):
        """List of SkyLineMember objects."""
        return self._members

    @property
    def polygons(self):
        """List of SkyLineMember polygons."""
        return [m.polygon for m in self.members]

    @classmethod
    def from_radec(cls, ra, dec, center=None, degrees=True):
        """See SphericalPolygon."""
        return SphericalPolygon.from_radec(ra, dec, center=center,
                                           degrees=degrees)

    @classmethod
    def from_cone(cls, ra, dec, radius, degrees=True, steps=16.0):
        """See SphericalPolygon."""
        return SphericalPolygon.from_cone(ra, dec, radius, degrees=degrees,
                                          steps=steps)

    @classmethod
    def from_wcs(cls, fitspath, steps=1, crval=None):
        """See SphericalPolygon."""
        return SphericalPolygon.from_wcs(fitspath, steps=steps, crval=crval)

    @classmethod
    def multi_union(cls, polygons, method='parallel'):
        """See SphericalPolygon."""
        return SphericalPolygon.multi_union(polygons, method=method)

    def union(self, other):
        """
        Updates *self* with the union of *self* and *other*.
        Skips *other* members that are already in *self*.

        Parameters
        ----------
        self: obj
            `SkyLine` instance to be updated.

        other: obj
            `SkyLine` instance to be added.

        Examples
        --------
        >>> s1 = SkyLine('image1.fits')
        >>> s2 = SkyLine('image2.fits')
        >>> s1.union(s2)  # s1 is updated

        See also
        --------
        sphere.polygon.SphericalPolygon.union
        
        """
        new_members = [m for m in other.members if m not in self.members]
        if len(new_members) == 0:
            return

        all_poly = self.multi_union(self.polygons +
                                    [m.polygon for m in new_members])

        self._update(all_poly, self.members + new_members)


# Overload parent class with following changes
#    a. add own attr
#    b. update self var in-place, no return

    def add_image(self,skyline):
        """Make composite SkyLine"""
        pass

    def compute_overlap(self, skyline):
        """Return sphere object with intersect of 2 skylines.

        Wrapper of sphere overlap method.

        """
        pass

    def find_intersection(self, skyline):
        """
        Return WCS object of overlap of 2 skylines.

        """
        pass

    def create_wcs(self):
        """Create WCS from SkyLine object.

        .. note:: Use stwcs to define a plane using multiple HSTWCS object

        Returns
        -------
        wcs: obj
            New HSTWCS objects.

        """
        pass

    def footprints(self):
        """Compute edges of skyline."""
        pass


def test():
    """Basic use case."""
    pass
