import math
import numpy
import pylab

# compute the positions of 2 WDs initially at rest as they
# infall into one-another, assuming point masses.

# we work in CGS units
G = 6.67428e-8        # cm^3 g^{-1} s^{-2}
M_sun = 1.98892e33    # g

# a simple class to serve as a container for the orbital information
# for the two stars
class stars:
    
    def __init__ (self, M_star1=1, M_star2=0.1,
                  x1_init=0.0, y1_init=0.0,
                  x2_init=0.0, y2_init=0.0):

        self.npts = -1
        self.maxpoints = 2000

        # star1 properties
        self.M_star1 = M_star1
        self.x1_init = x1_init
        self.y1_init = y1_init
        self.x_star1 = numpy.zeros(self.maxpoints)
        self.y_star1 = numpy.zeros(self.maxpoints)
        self.vx_star1 = numpy.zeros(self.maxpoints)
        self.vy_star1 = numpy.zeros(self.maxpoints)

        # star2 properties
        self.M_star2 = M_star2
        self.x2_init = x2_init
        self.y2_init = y2_init
        self.x_star2 = numpy.zeros(self.maxpoints)
        self.y_star2 = numpy.zeros(self.maxpoints)
        self.vx_star2 = numpy.zeros(self.maxpoints)
        self.vy_star2 = numpy.zeros(self.maxpoints)

        self.t = numpy.zeros(self.maxpoints)

    def integrate(self, dt, time):


        # allocate storage for R-K intermediate results
        # y[0:3] will hold the star1 info, y[4:7] will hold the star2 info
        k1 = numpy.zeros(8, numpy.float64)
        k2 = numpy.zeros(8, numpy.float64)
        k3 = numpy.zeros(8, numpy.float64)
        k4 = numpy.zeros(8, numpy.float64)

        y = numpy.zeros(8, numpy.float64)
        f = numpy.zeros(8, numpy.float64)

        t = 0.0

        # initial conditions

        # star 1
        y[0] = self.x1_init  # initial star1 x position
        y[1] = self.y1_init  # initial star1 y position

        y[2] = 0.0           # initial star1 x-velocity
        y[3] = 0.0           # initial star1 y-velocity

        # star 2
        y[4] = self.x2_init  # initial star2 x position
        y[5] = self.y2_init  # initial star2 y position

        y[6] = 0.0           # initial star2 x-velocity
        y[7] = 0.0           # initial star2 y-velocity

        self.x_star1[0] = y[0]
        self.y_star1[0] = y[1]

        self.vx_star1[0] = y[2]
        self.vy_star1[0] = y[3]

        self.x_star2[0] = y[4]
        self.y_star2[0] = y[5]

        self.vx_star2[0] = y[6]
        self.vy_star2[0] = y[7]

        self.t[0] = t

        n = 1
        while (n < self.maxpoints and t < time):

            f = self.rhs(t, y)
            k1[:] = dt*f[:]

            f = self.rhs(t+0.5*dt, y[:]+0.5*k1[:])
            k2[:] = dt*f[:]

            f = self.rhs(t+0.5*dt, y[:]+0.5*k2[:])
            k3[:] = dt*f[:]
            
            f = self.rhs(t+dt, y[:]+k3[:])
            k4[:] = dt*f[:]

            y[:] += (1.0/6.0)*(k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])

            t = t + dt

            self.x_star1[n] = y[0]
            self.y_star1[n] = y[1]

            self.vx_star1[n] = y[2]
            self.vy_star1[n] = y[3]

            self.x_star2[n] = y[4]
            self.y_star2[n] = y[5]

            self.vx_star2[n] = y[6]
            self.vy_star2[n] = y[7]

            self.t[n] = t

            n += 1


        self.npts = n



    def rhs(self,t,y):

        f = numpy.zeros(8, numpy.float64)

        # y[0] = x_star1, y[1] = y_star1, y[2] = vx_star1, y[3] = vy_star1
        # y[4] = x_star2, y[5] = y_star2, y[6] = vx_star2, y[7] = vy_star2

        # unpack
        x_star1 = y[0]
        y_star1 = y[1]

        vx_star1 = y[2]
        vy_star1 = y[3]

        x_star2 = y[4]
        y_star2 = y[5]

        vx_star2 = y[6]
        vy_star2 = y[7]


        # distance between stars
        r = numpy.sqrt((x_star2 - x_star1)**2 + (y_star2 - y_star1)**2)


        f[0] = vx_star1  # d(x_star1) / dt
        f[1] = vy_star1  # d(y_star1) / dt

        f[2] = -G*self.M_star2*(x_star1 - x_star2)/r**3  # d(vx_star1) / dt
        f[3] = -G*self.M_star2*(y_star1 - y_star2)/r**3  # d(vy_star1) / dt

        f[4] = vx_star2  # d(x_star2) / dt
        f[5] = vy_star2  # d(y_star2) / dt

        f[6] = -G*self.M_star1*(x_star2 - x_star1)/r**3  # d(vx_star2) / dt
        f[7] = -G*self.M_star1*(y_star2 - y_star1)/r**3  # d(vy_star2) / dt

        return f
    

def infall():

    # set the masses
    M_star1 = 1.19332959507310304e33
    M_star2 = 1.59123272493435018e33

    # set the semi-major axis of the stars
    a_star1 = 2.e10 - 3275397144.3388567
    a_star2 = 2.e10 + 2456352415.7779908

    # normalization
    x_center = 2.e10


    # create the star system object
    ss = stars(M_star1=M_star1, M_star2=M_star2,
               x1_init=a_star1, y1_init=0.0,
               x2_init=a_star2, y2_init=0.0)

    
    # set the timestep
    dt = 0.1        
    tmax = 33.0

    ss.integrate(dt,tmax)


    # output
    n = 0
    while (n < ss.npts):
        print ss.t[n], ss.x_star1[n]-x_center, ss.x_star2[n]-x_center
        n += 1

    
if __name__== "__main__":
    infall()


    
        
