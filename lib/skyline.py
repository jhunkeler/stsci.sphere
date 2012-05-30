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
from polygon import SphericalPolygon

__all__ = ['SkyLine']

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
        # Find SCI extensions
        with pyfits.open(fname) as pf:
            sci_extnum = tuple([i for i,ext in enumerate(pf) if
                                'SCI' in ext.name.upper()])
        assert len(sci_extnum) > 0, \
            'SkyLine cannot find SCI ext in {}.'.format(fname)

        # Initialize mosaic with first chip
        all_sph = SphericalPolygon.from_wcs(
            wcsutil.HSTWCS(fname, ext=sci_extnum[0]))

        # Mosaic all the chips
        for i in sci_extnum[1:]:
            prev_sph = all_sph.__copy__()
            
            cur_sph = SphericalPolygon.from_wcs(
                wcsutil.HSTWCS(fname, ext=i))
                
            all_sph = prev_sph.union(cur_sph)

        # SkyLine = final mosaic inheriting from SphericalPolygon
        SphericalPolygon.__init__(self, all_sph.points, all_sph.inside)

        # Add attribute to track images in SkyLine
        self._image_list = [fname]

    @property
    def image_list(self):
        """List of images in SkyLine."""
        return self._image_list

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
