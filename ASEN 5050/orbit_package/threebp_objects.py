import numpy as np

class Cr3bp():
    """ 
    This objet is used for work in the CR3BP, with built in funciotns for Jacobi integrals,
    potential functions and other utilities
    """
    def V_planar(self,x:float,y:float,mu:float)->float:
        """ 
        This function calculates the potential function value for a point in space
        for a given value of mu, assuming planar motion

        Args:
            x,y: Object position in 3d 
            mu: ratio of masses of larger bodies

        Returns:
            V: value of the potential funciton at that point
        """
        V = 0.5*(x**2+y**2) + (1-mu)/(np.sqrt((x+mu)**2+y**2)) + mu/np.sqrt((x-1+mu)**2+y**2)
        return V
    
    def V(self,x:float,y:float,z:float,mu:float)->float:
        """ 
        This function calculates the potential function value for a point in space
        for a given value of mu, assuming planar motion

        Args:
            x,y,z: Object position in 3d 
            mu: ratio of masses of larger bodies

        Returns:
            V: value of the potential funciton at that point
        """
        V = 0.5*(x**2+y**2) + (1-mu)/(np.sqrt((x+mu)**2+y**2+z**2)) + mu/np.sqrt((x-1+mu)**2+y**+z**2)
        return V
    
    def pVpx(self,x:float,y:float,z:float,mu:float)->float:
        """
        This function calculates the partial derivative of V
        with respect to x for a given position in 3D space

        """
        ppx = -mu*(mu+x-1)/((mu+x-1)**2+y**2+z**2)**(3/2) + ((mu-1)*(mu+x))/((mu+x)**2 + y**2+z**2)**(3/2) + x
        return ppx
    
    def pVpy(self,x:float,y:float,z:float,mu:float)->float:
        """
        This function calculates the partial derivative of V
        with respect to y for a given position in 3D space

        """
        ppy = -mu*y/((mu+x-1)**2+y**2+z**2)**(3/2) + ((mu-1)*y)/((mu+x)**2 + y**2+z**2)**(3/2) + y
        return ppy
    
    def pVpz(self,x:float,y:float,z:float,mu:float)->float:
        """
        This function calculates the partial derivative of V
        with respect to x for a given position in 3D space

        """
        ppz = -mu*z/((mu+x-1)**2+y**2+z**2)**(3/2) + ((mu-1)*z)/((mu+x)**2 + y**2+z**2)**(3/2) 
        return ppz
    
    def J_planar(self,x:float,y:float,vx:float,vy:float,mu:float)->float:
        """
        This function evaluates the Jacobi integral at a point in space with 
        a specified velocity value. This assumes purely planar motion and position

        Args:
            x,y: Position
            vx,vy: Velocity
            mu: mass ratio

        Returns:
            J: Integral Value
        """
        J = 0.5*(vx**2 + vy**2) - self.V_planar(x,y,mu)
        return J
