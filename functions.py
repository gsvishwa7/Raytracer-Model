# -*- coding: utf-8 -*-
"""
This script contains functions to observe or calculate physical properties of different optical elements

"""

import numpy as np
import raytracer as rt
import matplotlib.pyplot as plt

ray1 = rt.Ray(np.array([0.1,0,0]),np.array([0,0,1]))
ray2 = rt.Ray(np.array([-0.1,0,0]),np.array([0,0,1]))
rays = [ray1,ray2]

trial_ray1 = rt.Ray(p=np.array([20,0,0]),k=np.array([0,0,1]))
trial_ray2 = rt.Ray(p=np.array([-20,0,0]),k=np.array([0,0,1]))
trial_ray3 = rt.Ray(p=np.array([2,0,0]),k=np.array([0,0,1]))
trial_ray4 = rt.Ray(p=np.array([-2,0,0]),k=np.array([0,0,1]))
trial_rays = [trial_ray1, trial_ray2, trial_ray3, trial_ray4]

#New change 

def see_aberration():
    for trial_ray in trial_rays:
        for elem in rt.OpticalElement.elements:
            elem.propagate_ray(trial_ray)
        vertices = np.asarray(trial_ray.vertices())
        x = vertices[:,0]
        y = vertices[:,1]
        z = vertices[:,2]  
        
        if trial_ray == trial_ray1:
            #Output the marginal focus
            xx1 = vertices[1,0]
            point_z = vertices[1,:] + abs(xx1/trial_ray.k()[0])*(trial_ray.k())
            print('The marginal focus is at {}'.format(point_z))

        plt.figure(5)        
        plt.plot(z,x, c="b")
        plt.title('Spherical Aberration with 4 Rays')
        plt.xlabel('z')
        plt.ylabel('x')
    
    
def check_paraxial_focus():
    for ray in rays:
        for elem in rt.OpticalElement.elements:
            elem.propagate_ray(ray)
        vertices = np.asarray(ray.vertices())
        xx = vertices[-2,0]
        point1 = vertices[-2,:] + abs(xx/ray.k()[0])*(ray.k())
    print('The paraxial focus is at:{}'.format(point1))
    
    # Need more problems

def spherical_aberration_demo():
    demo_rays = []
    x_start = np.arange(0,22,2)
    for x in x_start:
        demo_rays.append(rt.Ray(p=np.array([x,0,0]),k=np.array([0,0,1])))
        demo_rays.append(rt.Ray(p=np.array([-x,0,0]),k=np.array([0,0,1])))
    for ray_d in demo_rays:
        for elem in rt.OpticalElement.elements:
            elem.propagate_ray(ray_d)
        vertices = np.asarray(ray_d.vertices())
        x = vertices[:,0]
        y = vertices[:,1]
        z = vertices[:,2]
        plt.figure(6)
        plt.plot(z,x, c="b")
        plt.title('Spherical Aberration')
        plt.xlabel('z')
        plt.ylabel('x')

def rms_calc():
    bundles=[]
    rms=[]
    values_j=[]
    number = np.arange(2,22,2)
    values_ds=[]
    for i in number:
        bundles.append(rt.RayBundle(radius=i, centre = np.array([0,0,0]), direction =np.array([0,0,1])))
        values_j.append(2*i)
        
    for bundle in bundles:
        rz=[]
        for ray in bundle.all_rays:
            for elem in rt.OpticalElement.elements:
                elem.propagate_ray(ray)
            vertices = np.asarray(ray.vertices())
            final_pos_x = vertices[-1,0]
            final_pos_y = vertices[-1,1]
            r_squared_val = (final_pos_x)**2 + (final_pos_y)**2
            rz.append(r_squared_val)
      
        rz = np.asarray(rz)
        rms.append(np.sqrt(np.sum(rz)/len(rz)))
    
    plt.figure(8)
    plt.plot(number,rms,c='g',label='RMS Spot Radius C2',marker='x')
     
    for j in values_j:
        ds=((588*10**(-6))*96.749)/(j)
        values_ds.append(ds)
        
    plt.plot(number,values_ds,c='k',linestyle='--',marker='.',label='Diffraction Limit')
    plt.legend()
    
plt.show()
