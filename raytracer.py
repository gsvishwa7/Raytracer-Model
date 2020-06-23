# -*- coding: utf-8 -*-
"""
Raytracer Module
By Girish S Vishwa
CID: 01218336

The values in all the scripts are set to replicate results obtained in the report
Some parts of some scripts are commented out as they are simply additional functions
These can be uncommented at the user's discretion 

"""

import numpy as np

class Ray:
    def __init__(self, p = np.array([0,0,0]), k= np.array([1,1,0])):
        if type(p) != np.ndarray or type(k) != np.ndarray:
            raise Exception('Enter a numpy array as input.')
        if  p.size != 3 or k.size != 3:
            raise Exception('Enter an array corresponding to a point/vector in 3D.')
        if np.linalg.norm(k) == 0:
            raise Exception('[0,0,0] cannot be a direction.')
        self.__vertices=[p]
        self.__p = p
        self.__k = k
        self.__intercepted = 0
    def p(self):
        return self.__p
    def k(self):
        return self.__k
    def append(self,new_p, new_k):
        if type(new_p) != np.ndarray or type(new_k) != np.ndarray:
            raise Exception('Enter a numpy array as input.')
        if  new_p.size != 3 or new_k.size != 3:
            raise Exception('Enter an array corresponding to a point/vector in 3D.')
        if np.linalg.norm(new_k) == 0:
            raise Exception('[0,0,0] cannot be a direction.')
        if not np.array_equal(new_p, self.__p):
            self.__vertices.append(new_p)
        self.__p = new_p
        self.__k = new_k
    def vertices(self):
        return self.__vertices
    def check_intercept(self):
        return self.__intercepted
    def set_intercept(self,new):
        self.__intercepted = new



class RayBundle:
    'Creates a bundle of rays with a step increment of 2mm for Radius'
    'This can be modified by tweaking the linspace conditions for r'
    def __init__(self, radius = 6, centre = np.array([-10,0,0]), direction = np.array([0,0,1])):    
        self.radius = radius
        self.direction = direction
        self.centre = centre
        self.all_rays=[]
        s = self.centre
        r = np.linspace(2,self.radius,(self.radius)/2)
        self.all_rays.append(Ray(p=self.centre,k=self.direction))
        for c in r:
            angles = np.arange(0,2*(np.pi),(np.pi/3)/c)
            for theta in angles:
                new_ray = Ray(p=np.array([c*np.cos(theta)+s[0],c*np.sin(theta)+s[1],0+s[2]]),k=self.direction)
                self.all_rays.append(new_ray)
            
    

class OpticalElement:
    elements = []
    i = 0 
    def __init__(self, z0, a0, n1, n2):
        self.z0 = z0
        self.a0 = a0
        self.n1 = n1
        self.n2 = n2
        self.update_element_order()
        
    def update_element_order(self):
        "Before propagting a ray, the ordering of the optical elements matters"
        "The Output Plane must always be the last object. This method ensures this."
        OpticalElement.elements.append(self)
        while OpticalElement.i < len(OpticalElement.elements)-1:
            if isinstance(OpticalElement.elements[OpticalElement.i], OutputPlane):
                OpticalElement.elements.append(OpticalElement.elements.pop(OpticalElement.i))
            else:
                OpticalElement.i += 1

    def propagate_ray(self, ray):
        "propagate a ray through the optical element"
        if isinstance(self, SphericalRefraction):
            self.intercept(ray)
            if ray.check_intercept() == 1:
                self.snell_refraction(ray)
        elif isinstance(self, PlaneRefraction):
            self.intercept(ray)
            if ray.check_intercept() == 1:
                #print('refracted')
                self.snell_refraction(ray)
        elif isinstance(self, SphericalMirror):
            self.intercept(ray)
            if ray.check_intercept() == 1:
                self.reflect(ray)
        elif isinstance(self, OutputPlane):
            self.intercept(ray)



class SphericalRefraction(OpticalElement):
    def __init__(self, z0=0, a0=0.5, n1 = 1, n2 = 1.5, R=25):
        if R == 0:
           raise Exception('Invalid Radius of Curvature')
        if a0 > abs(R):
            raise Exception('Invalid Aperture Radius')
        self.R = R
        OpticalElement.__init__(self, z0, a0, n1, n2)
        
    def intercept(self,ray):
        assert isinstance(ray, Ray), 'This object is not a ray.'
        #Define necessary variables first
        start_point = ray.p()
        mag_k = np.linalg.norm(ray.k())
        k_hat = ray.k()/mag_k
        r = start_point - np.array([0,0,self.R+self.z0])
        discriminant = (np.dot(r,k_hat))**2-((np.linalg.norm(r))**2-(self.R)**2)
        x = start_point[0]
        y = start_point[1]
        z = start_point[2]
        if x**2 + y**2 + (z-(self.R + self.z0))**2 > (self.R)**2:
            if np.dot(r,k_hat) >= 0:
                ray.set_intercept(0)
                return None
            else:
                 if discriminant < 0:
                     ray.set_intercept(0)
                     return None
                 else:
                     if self.R>0:
                         k2 = -np.dot(r,k_hat)-np.sqrt(discriminant)
                         interception_point = start_point + k2*k_hat
                         ray.append(interception_point,ray.k())
                     elif self.R<0:
                         k2 = -np.dot(r,k_hat)+np.sqrt(discriminant)
                         interception_point = start_point + k2*k_hat
                         ray.append(interception_point,ray.k())
                     ray.set_intercept(1)
                     return interception_point
        elif x**2 + y**2 + (z-(self.R + self.z0))**2 < (self.R)**2:
            if self.R > 0:
                raise Exception('Interception does not correspond to +ve R')
            k2 = -np.dot(r,k_hat)+np.sqrt(discriminant)
            interception_point = start_point + k2*k_hat
            ray.append(interception_point,ray.k())
            ray.set_intercept(1)
            return interception_point
      

    def snell_refraction(self,ray):
        normal_unit = (np.array([0,0,self.z0+self.R])-ray.p())/self.R
        incident_direction = ray.k()        
        #split the incident ray direction into the two relevant components
        ray_parallel = (np.dot(normal_unit, incident_direction))*normal_unit
        ray_perpendicular = incident_direction - ray_parallel
        theta1 = np.arccos((np.dot(incident_direction, normal_unit))/np.linalg.norm(incident_direction))
        if np.sin(theta1) > self.n2/self.n1:
            return None
        elif theta1 == 0 or theta1 == np.pi:
            new_direction = incident_direction
            ray.append(ray.p(),new_direction)
        else:
            theta2 = np.arcsin((self.n1/self.n2)*np.sin(theta1))
            new_parallel_scale = np.linalg.norm(ray_perpendicular)/(np.tan(theta2))
            new_direction = ray_perpendicular + ((new_parallel_scale)/(np.linalg.norm(ray_parallel)))*ray_parallel
            ray.append(ray.p(),new_direction)
        
        
            
class PlaneRefraction(OpticalElement):
    def __init__(self, z0=0, a0=1000, n1 = 1, n2 = 1.5, b0 = 1000):
        #a0 sets the extent of the plane along the x-axis 
        #b0 sets the extent of the plane along the y-axis
        self.b0 = b0
        OpticalElement.__init__(self, z0, a0, n1, n2)
    def intercept(self, ray):
        ray.set_intercept(0)
        if self.z0 - ray.p()[2] >0 and ray.k()[2] < 0 or self.z0 - ray.p()[2] <0 and ray.k()[2] >0 or ray.k()[2] == 0:
            ray.set_intercept(0)
            return None
        else:
            new_coord = (abs((self.z0 - ray.p()[2])/(ray.k()[2])))*ray.k() + ray.p()
            if -self.b0 < new_coord[1] < self.b0 and -self.a0 < new_coord[0] < self.a0:
                ray.append(new_coord,ray.k())
                ray.set_intercept(1)
            else:
                ray.set_intercept(0)
                return None
            
    def snell_refraction(self, ray):
        incident_direction = ray.k()
        if incident_direction[2] > 0:
            normal_unit = np.array([0,0,1])
        else:
            normal_unit = np.array([0,0,-1])
        ray_parallel = (np.dot(normal_unit, incident_direction))*normal_unit
        ray_perpendicular = incident_direction - ray_parallel 
        theta1 = np.arccos((np.dot(incident_direction, normal_unit))/(np.linalg.norm(incident_direction)))
        if np.sin(theta1) > self.n2/self.n1:
            return None
        elif theta1 == 0 or theta1 == np.pi:
            new_direction = incident_direction
            ray.append(ray.p(),new_direction)
        else:
            theta2 = np.arcsin((self.n1/self.n2)*np.sin(theta1))
            new_parallel_scale = np.linalg.norm(ray_perpendicular)/(np.tan(theta2))
            new_direction = ray_perpendicular + ((new_parallel_scale)/(np.linalg.norm(ray_parallel)))*ray_parallel
            ray.append(ray.p(),new_direction)

        
        
        
class OutputPlane(OpticalElement):
    def __init__(self,z0):
        self.z0=z0
        OpticalElement.__init__(self,z0,a0=100,n1=1,n2=1)
    def intercept(self, ray):
        if self.z0 - ray.p()[2] >0 and ray.k()[2] < 0 or self.z0 - ray.p()[2] <0 and ray.k()[2] >0 or ray.k()[2] == 0:
            return None
        else:
            new_coord = (abs((self.z0 - ray.p()[2])/(ray.k()[2])))*ray.k() + ray.p()
            ray.append(new_coord,ray.k())
            

class SphericalMirror(OpticalElement):
    def __init__(self, z0=0, a0=0.5, n1 = 1, n2 = 1.5, R=25):
        if R == 0:
           raise Exception('Invalid Radius of Curvature')
        if a0 > abs(R):
            raise Exception('Invalid Aperture Radius')
        self.R = R
        OpticalElement.__init__(self, z0, a0, n1, n2)
        
    def intercept(self,ray):
        SphericalRefraction.intercept(self,ray)

    def reflect(self,ray):
        "this method simply involves flipping the component of the incoming ray which is parallel to the normal"
        normal_unit = (np.array([0,0,self.z0+self.R])-ray.p())/self.R
        #split the incident ray direction into the two relevant components
        incident_direction = ray.k()
        ray_parallel = (np.dot(normal_unit, incident_direction))*normal_unit
        ray_perpendicular = incident_direction - ray_parallel
        new_direction = ray_perpendicular - ray_parallel
        ray.append(ray.p(),new_direction)
        
        
    
        
        
