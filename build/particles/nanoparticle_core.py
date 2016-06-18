""":module:: nanoparticle_core
   :platform: Unix, Windows
   :synopsis: checks if a point is within the unit cell of a bcc lattice

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February 2016
"""

from __future__ import division

import numpy as np
import numpy.linalg as la

class NanoBcc(object):
    """Defines points that fall within the unit cell of a bcc lattice
    """

    def __init__(self, radius):
        """The constructor
        
        particles are assumed to have diameter 1. All units are relative to 
        this one
        
        :param radius: radius of unit cell
        """
        
        cfac = 2/np.sqrt(3)
        
        # those are the 14 nearest neighbors of the (0,0,0) point
        self.pnt = cfac*np.array([
        [0.5, 0.5, 0.5],
        [0.5, 0.5, -0.5],
        [0.5, -0.5, 0.5],
        [-0.5, 0.5, 0.5],
        [0.5, -0.5, -0.5],
        [-0.5, 0.5, -0.5],
        [-0.5, -0.5, 0.5],
        [-0.5, -0.5, -0.5],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0]]
        )
        
        # we now compute the coordinates of the planes
        # the first three numbers are point coordinates, fourth is intercept
        self.planes = np.zeros([self.pnt.shape[0], 4])
        
        for ind, bar in enumerate(self.pnt):
            cnorm = la.norm(bar)
            # compute normals
            self.planes[ind][:-1] = bar[:]/cnorm     
            # compute intercept
            self.planes[ind][-1] = - cnorm*(radius-0.5)
            
    def check_point(self, pnt):
        """Checks whether a  point is within the unit cell
        
        :param pnt: given point
        :return: True if the point is within the unit cell
        :rtype: bool
        """
        
        for ind in range(self.planes.shape[0]):
            if np.dot(self.planes[ind,:-1], pnt) + self.planes[ind, -1] > 0:
                return False
        return True


class NanoFcc(object):
    """Defines points that fall within the unit cell of a fcc lattice
    """

    def __init__(self, radius):
        """The constructor
        
        particles are assumed to have diameter 1. All units are relative to 
        this one
        
        :param radius: radius of unit cell
        """
        
        cfac = np.sqrt(2)/2
        
        # those are the 18 nearest neighbors of the (0,0,0) point
        self.pnt = cfac*np.array([
        [0.5, 0.5, 0.0],
        [-0.5, -0.5, 0.0],
        [0.5, -0.5, 0.0],
        [-0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [-0.5, 0.0, -0.5],
        [-0.5, 0.0, 0.5],
        [0.5, 0.0, -0.5],
        [0.0, 0.5, 0.5],
        [0.0, -0.5, -0.5],
        [0.0, 0.5, -0.5],
        [0.0, -0.5, 0.5],
        [1.0, 0.0 ,0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0]]
        )
        
        # we now compute the coordinates of the planes
        # the first three numbers are point coordinates, fourth is intercept
        self.planes = np.zeros([self.pnt.shape[0], 4])
        
        for ind, bar in enumerate(self.pnt):
            cnorm = la.norm(bar)
            # compute normals
            self.planes[ind][:-1] = bar[:]/cnorm     
            # compute intercept
            self.planes[ind][-1] = - cnorm*(radius-0.5)
            
    def check_point(self, pnt):
        """Checks whether a  point is within the unit cell
        
        :param pnt: given point
        :return: True if the point is within the unit cell
        :rtype: bool
        """
        
        for ind in range(self.planes.shape[0]):
            if np.dot(self.planes[ind,:-1], pnt) + self.planes[ind, -1] > 0:
                return False
        return True



if __name__ == '__main__':
    
    print('Unit Test for NanoBcc')
    
    uc_size = 4
    
    print('We compute a unit cell of size %d', uc_size)
    
    uc_bcc = NanoBcc(uc_size)

    print('(0, 0, 0), (3, 0, 0), (0, 3, 0), (0, 3, 1)) are within unit cell')

    cfac = 2.0/np.sqrt(3)

    vec0 = np.array([0.0, 0.0, 0.0])
    vec1 = cfac*np.array([3.0, 0.0, 0.0]) 
    vec2 = cfac*np.array([0.0, 3.0, 0.0])
    vec3 = cfac*np.array([0.0, 0.0, 3.0])
    
    vec_all = [vec0, vec1, vec2, vec3]
    
    for vec in vec_all:
        if not uc_bcc.check_point(vec):
            print(vec)
            raise ValueError('unit test did not pass')
    
    
    print('(4, 0, 0), (0, 4, 0), (0, 0, 4)) are outside unit cell')

    vec1 = cfac*np.array([4.0, 0.0, 0.0]) 
    vec2 = cfac*np.array([0.0, 4.0, 0.0])
    vec3 = cfac*np.array([0.0, 0.0, 4.0])   
    
    vec_all = [vec1, vec2, vec3]    
    
    for vec in vec_all:
        if uc_bcc.check_point(vec):
            print(vec)
            raise ValueError('unit test did not pass')
            
    print('UNIT TEST PASSED')
