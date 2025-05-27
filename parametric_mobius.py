import numpy as np
import matplotlib.pyplot as plt


class MobiusStrip:
    def __init__(self, R, w, n):
        self.R = R
        self.w = w
        self.n = n

        """
        the above are the Mobius strip parameters initialized.
        R: Radius from center to the middle of the strip or is radius of base circle
        w: Width of the strip
        n: Resolution (number of points in u and v)
        """

        self.u = np.linspace(0, 2*np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)

        """
        the above are the surface parameters
        u is a parameter representing the angle around the base circle of the strip.
        v is a parameter representing the distance along the strip
        """
        
        self.U, self.V = np.meshgrid(self.u, self.v)

    
        #the above U,V are the coordinate grids of u,v data
        

        self.X, self.Y, self.Z = self.generate_xyz()
    
    def generate_xyz(self):
        #to generate x,y,z coordinte grid or 3D mesh using parametric mesh grid
        
        R = self.R
        U = self.U
        V = self.V

        X = (R + (V * np.cos(U/2))) * np.cos(U)
        Y = (R + (V * np.cos(U/2))) * np.sin(U)
        Z = V * np.sin(U/2)
        
        # the above are the parametric equations of the mobius strip

        return (X,Y,Z)
    
    def generate_mobius_plot(self):
        #using matplotlib to plot 3D visualization of mobius strip

        ax = plt.axes(projection="3d")

        ax.plot_surface(self.X, self.Y, self.Z,color='magenta', edgecolor='k')
        #edge color, the color of the turn in the strip, k for black

        ax.set_title("Mobius Strip")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.show()
    
    def calculate_surface_area_of_strip(self):
        """
        1. Approximated by integrating(summing up) over magnitude of the cross product 
        of first order partial derivates of small(differential) parametric surface w.r.t u,v.
        2. The first order partial derivatives represent the tangents to differential strip surface
        3. Because the cross product of tangents give normal vector to surface
        """

        du = self.u[1] - self.u[0]
        dv = self.v[1] - self.v[0]

        #calculating partial derivatives
        dX_du, dX_dv = np.gradient(self.X, du, dv, edge_order=2)
        dY_du, dY_dv = np.gradient(self.Y, du, dv, edge_order=2)
        dZ_du, dZ_dv = np.gradient(self.Z, du, dv, edge_order=2)

        #cross product to get normal_vectors of all differential surfaces

        normal_vectors = np.cross(
            np.stack((dX_du, dY_du, dZ_du), axis=-1),
            np.stack((dX_dv, dY_dv, dZ_dv), axis=-1)
        )

        #taking the magnitude of the product

        magnitudes = np.linalg.norm(normal_vectors, axis=-1)
        
        #calculating surface area by summing up

        area = np.sum(magnitudes) * du * dv
        return area

    def calculate_edge_length(self):
        #by summing up small edge lengths of all differential surfaces

        u_vals = self.u
        v_edge = self.w/2
        R = self.R
    
        #again using parametric surface
        x_edge = (R + v_edge * np.cos(u_vals/2)) * np.cos(u_vals)
        y_edge = (R + v_edge * np.cos(u_vals/2)) * np.sin(u_vals)
        z_edge = v_edge * np.sin(u_vals/2)

        dx = np.diff(x_edge)
        dy = np.diff(y_edge)
        dz = np.diff(z_edge)
        segment_lengths = np.sqrt(dx**2 + dy**2 + dz**2)
        length = np.sum(segment_lengths)
        return length



mobius = MobiusStrip(R=5, w=2, n=200)


surface_area = mobius.calculate_surface_area_of_strip()
edge_length = mobius.calculate_edge_length()

print("Surface Area of mobius strip given: {}".format(surface_area))
print("Edge Length of the mobius strip: {}".format(edge_length))

mobius.generate_mobius_plot()



        
