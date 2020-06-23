# -*- coding: utf-8 -*-
"""
This script can be used to observe the behaviour of a Spherical Mirror

"""
import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt
import functions as f

rt.OpticalElement.elements = []
output = rt.OutputPlane(90)
o = rt.SphericalMirror(z0=120, R=-50)
B = rt.RayBundle(radius=10, centre = np.array([0,0,0]), direction = np.array([0,0,1]))

for ray2 in B.all_rays:
    for elem in rt.OpticalElement.elements:
        elem.propagate_ray(ray2)
        
    vertices2 = np.asarray(ray2.vertices())
    y2 = vertices2[:,1]
    x2 = vertices2[:,0]
    z2 = vertices2[:,2]
    plt.figure(4)
    plt.plot(z2,x2, c="b")
    plt.title('Spherical Surface Refraction')
    plt.xlabel('z')
    plt.ylabel('x')

f.spherical_aberration_demo()
f.check_paraxial_focus()
#f.see_aberration()
