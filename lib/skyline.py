"""
skyline -- Manage outlines on the sky

This module provides support for working with footprints on the sky.
Primary use case would use the following generalized steps::

1. Initialize SkyLine objects for each input image. This object would be the 
    union of all the input image's individual chips WCS footprints.
2. Determine overlap between all images. The determination would employ a 
    recursive operation to return the extended list of all overlap values 
    computed as [img1 vs [img2,img3,...,imgN],img2 vs [img3,...,imgN],...]
3. Select the pair with the largest overlap, or the pair which produces the 
    largest overlap with the first input image. This defines the initial 
    reference SkyLine object.
4. Perform some operation on the 2 images: for example, match sky in intersecting
    regions, or aligning second image with the first (reference) image.
5. Update the second image, either apply the sky value or correct the WCS, then 
    generate a new SkyLine object for that image.
6. Create a new reference SkyLine object as the union of the initial reference
    object and the newly updated SkyLine object.
7. Repeat Steps 2-6 for all remaining input images.

This process will work reasonably fast as most operations are performed using
the SkyLine objects and WCS information solely, not image data itself.
"""
import pyfits

import sphere

class SkyLine(object):
    def __init__(self,fname):
        pass
        
    def add_image(self,skyline):
        pass
    
    def compute_overlap(self, skyline):
        pass
    
    def find_intersection(self, skyline):
        pass
    
    def within(self, pos):
        pass
    
    def create_wcs(self):
        pass
        
