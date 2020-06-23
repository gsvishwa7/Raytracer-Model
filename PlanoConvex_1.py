# -*- coding: utf-8 -*-
"""
This script is to study the properties of a plano-convex lens with the planar side facing the collimated light beam
"""
import numpy as np
import matplotlib.pyplot as plt
import raytracer as rt
import functions as f

rt.OpticalElement.elements = []
output = rt.OutputPlane(201.7488)
p1 = rt.PlaneRefraction(a0=10000000,b0=10000000,z0=100,n1=1,n2=1.5168)
o = rt.SphericalRefraction(z0=105, R=-50, n1=1.5168,n2=1)
C = rt.RayBundle(radius=10, centre = np.array([0,0,0]), direction =np.array([0,0,1]))

r = []

for ray in C.all_rays:

    for elem in rt.OpticalElement.elements:
        elem.propagate_ray(ray)
    vertices = np.asarray(ray.vertices())
    x = vertices[:,0]
    y = vertices[:,1]
    z = vertices[:,2]
    init_pos_x = vertices[0,0]
    init_pos_y = vertices[0,1]
    final_pos_x = vertices[-1,0]
    final_pos_y = vertices[-1,1]
    r_squared_val = (final_pos_x)**2 + (final_pos_y)**2
    r.append(r_squared_val)
    
    plt.figure(1)
    plt.plot(z,x, c="b")
    plt.ylim((-20,20))
    plt.title('Plano-Convex Lens with Planar Side Facing Collimated Light Beam')
    plt.xlabel('z')
    plt.ylabel('x')
        
    plt.figure(2)
    plt.scatter(init_pos_x,init_pos_y, c="b")
    plt.axis('equal') 
    plt.title('x-y plane at z=0')
    plt.xlabel('x')
    plt.ylabel('y')
    
    plt.figure(3)
    plt.scatter(final_pos_x, final_pos_y, c = "b")
    plt.axis('equal')
    plt.title('x-y plane at paraxial focus')
    plt.xlabel('x')
    plt.ylabel('y')

##The RMS value will be correct only if the output plane is placed at the paraxial focus
r = np.asarray(r)
rms = np.sqrt(np.sum(r)/len(r))
print('The RMS spot radius is {}'.format(rms))

f.check_paraxial_focus()
f.rms_calc()
#f.spherical_aberration_demo()


plt.show()