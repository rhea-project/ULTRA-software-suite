## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil-CST
# create_CST_from_airfoil.py
# 
# Created:  Jun 2020, S.Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from SUAVE.Core import Data
from SUAVE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil import import_airfoil_Lednicer as airfoil_geo
import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import minimize
from math           import pi, cos, factorial


# ------------------------------------------------------------
#  Creates CST coefficient for a given airfoil
# ------------------------------------------------------------

## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil-CST
class CST_from_aifroil:

    def __init__(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
            Possible shape coefficients:
                                        N1 = 0.5    N2 = 1.0   - NACA type round nose and pointed aft end airfoil
                                        N1 = 0.5    N2 = 0.5   - elliptic airfoil
                                        N1 = 1.0    N2 = 1.0   - biconvex airfoil

        Source:
        None
        
        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """

        self.shape_coefficients    = Data()
        self.shape_coefficients.N1 = 0.5
        self.shape_coefficients.N2 = 1.0
        self.shape_coefficients.dz = 0.0                                          # trailing edge gap from the chord line


    def split_airfoil(self,x1,y1,N):
        """Splits airfoil point into upper and lower surfaces
        
        Assumptions:
        Airfoil file in Lednicer format

        Source:
        None                 

        Inputs:
        x1 - airfoil X-coordinates
        N  - number of airfoil points

        Outputs:
        xl,xu -x-coordinates of upper and lower surfaces

        Properties Used:
        N/A
        """      
        center_loc = np.where(x1 == 0)  # Used to separate upper and lower surfaces
        center_loc = center_loc[0][0]

        xl = np.zeros(center_loc+1)
        xu = np.zeros(N-center_loc)
        yl = np.zeros(center_loc+1)
        yu = np.zeros(N-center_loc)

        for i in range(len(xl)):
            xl[i] = x1[i]            # Lower surface x-coordinates
            yl[i] = y1[i]
        for i in range(len(xu)):
            xu[i] = x1[i + center_loc]   # Upper surface x-coordinates
            yu[i] = y1[i + center_loc]
        return xl, xu, yl,yu


    def create_CST_from_airfoil(self,filename,order,dz):
        """
        Runs an algorythm to compute CST coefficients for a given airfoil
                    
        Assumptions:
        Airfoil file in Lednicer format

        Source:
            B. Kulfan CST - Universal Parametric Geometry Representation Method With Applications to Supersonic Aircraft
            Fourth International Conference on Flow Dynamics Sendai International Center Sendai, Japan September 26-28, 2007

        Inputs:
        filename   <string>

        Outputs:
        data       numpy array with airfoil data

        Properties Used:
        N/A
        """
        
        # unpack inputs
        N1 = self.shape_coefficients.N1
        N2 = self.shape_coefficients.N2        

        # import airfoil points       
        x0, y0, N = airfoil_geo(filename)
        
        # Split airfoil coordinates into upper and lower surfaces
        xll, xuu, yll, yuu = self.split_airfoil(x0,y0,N)
        N = int(N/2)

        # Minimize the lower airfoil surface
        def min_error(x):

            wl = []
            for i in range(order):
                wl.append(x[i])

            yl = self.create_CST(N1,N2,wl,xll,-dz)

            f = 0
            for i in range(N-1):
                f += abs(yl[i] - yll[i])**2            
            return f
        x = np.ones(order)

        res1 = minimize(min_error, x, method ='SLSQP', options={'ftol': 1e-9, 'disp': False})

        # Minimize the upper airfoil surface
        def min_error(x):

            wu = []
            for i in range(order):
                wu.append(x[i])

            yu = self.create_CST(N1,N2,wu,xuu,dz)

            f = 0
            for i in range(N-1):
                f += abs(yu[i] - yuu[i])**2            
            return f
        x = np.ones(order)

        res2 = minimize(min_error, x, method ='SLSQP', options={'ftol': 1e-9, 'disp': False})

        # Combine results for both surfaces into one       
        res = np.concatenate((res1.x,res2.x), axis=0)
        
        return res       
 
    def create_CST(self,N1,N2,w,x,dz):
        """
        CST generation algorythm
                    
        Assumptions:
        Airfoil file in Lednicer format

        Source:
            B. Kulfan CST - Universal Parametric Geometry Representation Method With Applications to Supersonic Aircraft
            Fourth International Conference on Flow Dynamics Sendai International Center Sendai, Japan September 26-28, 2007

        Inputs:
        filename   <string>

        Outputs:
        data       numpy array with airfoil data

        Properties Used:
        N/A
        """

        # Class function; taking input of N1 and N2
        C = np.zeros(len(x))
        for i in range(len(x)):
            C[i] = x[i]**N1*((1-x[i])**N2)

        # Shape function; using Bernstein Polynomials
        n = len(w) - 1  # Order of Bernstein polynomials

        K = np.zeros(n+1)
        for i in range(0, n+1):
            K[i] = factorial(n)/(factorial(i)*factorial(n-i))

        S = np.zeros(len(x))
        for i in range(len(x)):
            S[i] = 0
            for j in range(0, n+1):
                S[i] += w[j]*K[j]*x[i]**(j) * (1-x[i])**(n-j)

        # Calculate y output
        y = np.zeros(len(x))
        for i in range(len(y)):
            y[i] = C[i] * S[i] + x[i] * dz

        return y
