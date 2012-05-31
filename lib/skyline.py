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

        self._members = poly_list

        # Put mosaic of all the chips in SkyLine
        mosaic = SphericalPolygon.multi_union(self.polygons)
        self._points = mosaic.points
        self._inside = mosaic.inside

    @property
    def members(self):
        """List of SkyLineMember objects."""
        return self._members

    @property
    def polygons(self):
        """List of SkyLineMember polygons."""
        return [m.polygon for m in self.members]

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

    def within(self, pos):
        """
        Return bool if pos is in skyline or not.

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

    # Allow disjoint?

    pass
