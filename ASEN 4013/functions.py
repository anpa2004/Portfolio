import numpy as np
from pypdf import PdfReader
import re
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from collections import OrderedDict

def remove_letters_re(input_string):
    # The pattern '[a-zA-Z]' matches any uppercase or lowercase letter.
    # re.sub() replaces all occurrences of the pattern with an empty string.
    str1 = re.sub(r'[a-zA-Z]','',input_string)
    str2 = re.sub(r'[,/`()_¼:]', '', str1)
    return str2

def extract_appl(SI = True):
    # Create a PdfReader object by providing the path to your PDF file
    reader = PdfReader('SM 06_Mattingly_AppL.pdf')

    # Get the number of pages
    num_pages = len(reader.pages)

    if SI:
        TK = []
        TC = []
        hf0 = []
        Prf0 = []
        hf169 = []
        Prf169 = []
        hf338 = []
        Prf338 = []
        hf507 = []
        Prf507 = []
        hf676 = []
        Prf676 = []

        for i in range(num_pages):
            if i>0 and i<=10:
                page = reader.pages[i]
                text = page.extract_text()
                num_list = remove_letters_re(text)
                splat = num_list.split('\n')
                for x in splat:
                    if len(x)>55:
                        y = x.split(' ')
                        TK.append(float(y[0]))
                        TC.append(float(y[1]))
                        hf0.append(float(y[2]))
                        Prf0.append(float(y[3]))
                        hf169.append(float(y[4]))
                        Prf169.append(float(y[5]))
                        hf338.append(float(y[6]))
                        Prf338.append(float(y[7]))
                        hf507.append(float(y[8]))
                        Prf507.append(float(y[9]))
                        hf676.append(float(y[10]))
                        Prf676.append(float(y[11]))     
        return list(OrderedDict.fromkeys(TK)),list(OrderedDict.fromkeys(TC)), list(OrderedDict.fromkeys(hf0)),list(OrderedDict.fromkeys(Prf0)),list(OrderedDict.fromkeys(hf169)),list(OrderedDict.fromkeys(Prf169)),list(OrderedDict.fromkeys(hf338)),list(OrderedDict.fromkeys(Prf338)),list(OrderedDict.fromkeys(hf507)),list(OrderedDict.fromkeys(Prf507)),list(OrderedDict.fromkeys(hf676)),list(OrderedDict.fromkeys(Prf676))


    else:
        TR = []
        TF = []
        hf0 = []
        Prf0 = []
        hf169 = []
        Prf169 = []
        hf338 = []
        Prf338 = []
        hf507 = []
        Prf507 = []
        hf676 = []
        Prf676 = []

        for i in range(num_pages):
            if i>10:
                page = reader.pages[i]
                text = page.extract_text()
                num_list = remove_letters_re(text)
                splat = num_list.split('\n')
                for x in splat:
                    if len(x)>55:
                        y = x.split(' ')
                        TR.append(float(y[0]))
                        TF.append(float(y[1]))
                        hf0.append(float(y[2]))
                        Prf0.append(float(y[3]))
                        hf169.append(float(y[4]))
                        Prf169.append(float(y[5]))
                        hf338.append(float(y[6]))
                        Prf338.append(float(y[7]))
                        hf507.append(float(y[8]))
                        Prf507.append(float(y[9]))
                        hf676.append(float(y[10]))
                        Prf676.append(float(y[11]))     
        return list(OrderedDict.fromkeys(TR)),list(OrderedDict.fromkeys(TF)), list(OrderedDict.fromkeys(hf0)),list(OrderedDict.fromkeys(Prf0)),list(OrderedDict.fromkeys(hf169)),list(OrderedDict.fromkeys(Prf169)),list(OrderedDict.fromkeys(hf338)),list(OrderedDict.fromkeys(Prf338)),list(OrderedDict.fromkeys(hf507)),list(OrderedDict.fromkeys(Prf507)),list(OrderedDict.fromkeys(hf676)),list(OrderedDict.fromkeys(Prf676))

def linear_interp(x1,y1,x2,y2,xq):
    """  
    Linearly interpolate between two points x1,x2 at xq
    """
    m = (y2- y1)/(x2-x1)
    yq = y2 + m*(xq-x2)
    return yq

def cubic_spline_interpolation(x,y,xq):
    """ 
    This function does interpolation for one or several query points xq

    Args:
        x: list of x data
        y: list of y data
        xq: point or list of points to interpolate at
    
    Returns:
        yq: single value or list of values evaluated at the query point using cubic spline interpolation
    """

    cs = CubicSpline(x,y,bc_type='natural')
    yq = cs(xq)
    return yq

def find_f(Tt3,Tt4,etab,hpr,tol=1e-6,nmax=100,SI = True,K = True,pr = False,R=True):
    """ 
    Interpolation and Iteration scheme for fuel air ratio and enthalpy values

    Args:
        Tt3: Total Temperature at station 3
        Tt4: Total temperature at station 4
        etab: burner efficiency [0,1)
        hpr: Fuel specific heat thing
        tol: Iteration stopping tolerance, default to 1e-6
        nmax: Max iteration number
        SI: If to use SI units (bool). Must match Tt3 and Tt4. Defaults to True
        K: If units of K for Tt3, Tt4 and subsequent calcs, Defaults to True (False for Celsius)
        R: If units of R for Tt3, Tt4 and subsequent calcs, Defaults to True (False for Ferenheit)
        pr: If to print iteration process, defaults to False
    
    Returns:
        f: Fuel air ratio
        ht4: Enthalpy at station 4
    """
    flist = [0,0.0169,0.0338,0.0507,0.0676]
    if SI:
        TK,TC, hf0,Prf0,hf169,Prf169,hf338,Prf338,hf507,Prf507,hf676,Prf676 = extract_appl(SI=True)
        if K:
            # Initialization
            ht3 = cubic_spline_interpolation(TK,hf0,Tt3)
            f0 = flist[0]

            # Initial guess of f
            fn1 = flist[1] 

            # Interpolated values of h at Tt4 for each value of f
            hlist = [cubic_spline_interpolation(TK,hf0,Tt4),cubic_spline_interpolation(TK,hf169,Tt4),cubic_spline_interpolation(TK,hf338,Tt4),cubic_spline_interpolation(TK,hf507,Tt4),cubic_spline_interpolation(TK,hf676,Tt4)]
            
            # Initial guess of ht4 given initial guess of f
            ht4 = hlist[1]
            fn = (ht4-ht3)/(etab*hpr - ht4)
            n = 0

            # Iterating until the stopping conditions are met
            while np.abs(fn-fn1)>tol and n<nmax:
                if pr:
                    print(f"f = {fn}, ht4 = {ht4}, n = {n}")

                # Set new f_{n-1}
                fn1 = fn
                # Increase n to prevent infinite looping
                n+=1

                # Calculate ht4 at new value of f_{n-1}
                ht4 = cubic_spline_interpolation(flist,hlist,fn1)

                # Calculate new f_n given new ht4 value
                fn = (ht4 - ht3)/(etab*hpr - ht4)

                # Stopping infinite loop
                if n == nmax:
                    print('WARNING- Max iterations met')
            return fn,ht4
        else:
            ht3 = cubic_spline_interpolation(TC,hf0,Tt3) 
            f0 = flist[0]

            # Initial guess of f
            fn1 = flist[1] 

            # Interpolated values of h at Tt4 for each value of f
            hlist = [cubic_spline_interpolation(TC,hf0,Tt4),cubic_spline_interpolation(TC,hf169,Tt4),cubic_spline_interpolation(TC,hf338,Tt4),cubic_spline_interpolation(TC,hf507,Tt4),cubic_spline_interpolation(TC,hf676,Tt4)]
            
            # Initial guess of ht4 given initial guess of f
            ht4 = hlist[1]
            fn = (ht4-ht3)/(etab*hpr - ht4)
            n = 0

            # Iterating until the stopping conditions are met
            while np.abs(fn-fn1)>tol and n<nmax:
                if pr:
                    print(f"f = {fn}, ht4 = {ht4}, n = {n}")

                # Set new f_{n-1}
                fn1 = fn
                # Increase n to prevent infinite looping
                n+=1

                # Calculate ht4 at new value of f_{n-1}
                ht4 = cubic_spline_interpolation(flist,hlist,fn1)

                # Calculate new f_n given new ht4 value
                fn = (ht4 - ht3)/(etab*hpr - ht4)

                # Stopping infinite loop
                if n == nmax:
                    print('WARNING- Max iterations met')
            return fn,ht4
    else:
        TK,TC, hf0,Prf0,hf169,Prf169,hf338,Prf338,hf507,Prf507,hf676,Prf676 = extract_appl(SI=False)
        if R:
            # Initialization
            ht3 = cubic_spline_interpolation(TK,hf0,Tt3)
            f0 = flist[0]

            # Initial guess of f
            fn1 = flist[1] 

            # Interpolated values of h at Tt4 for each value of f
            hlist = [cubic_spline_interpolation(TK,hf0,Tt4),cubic_spline_interpolation(TK,hf169,Tt4),cubic_spline_interpolation(TK,hf338,Tt4),cubic_spline_interpolation(TK,hf507,Tt4),cubic_spline_interpolation(TK,hf676,Tt4)]
            
            # Initial guess of ht4 given initial guess of f
            ht4 = hlist[1]
            fn = (ht4-ht3)/(etab*hpr - ht4)
            n = 0

            # Iterating until the stopping conditions are met
            while np.abs(fn-fn1)>tol and n<nmax:
                if pr:
                    print(f"f = {fn}, ht4 = {ht4}, n = {n}")

                # Set new f_{n-1}
                fn1 = fn
                # Increase n to prevent infinite looping
                n+=1

                # Calculate ht4 at new value of f_{n-1}
                ht4 = cubic_spline_interpolation(flist,hlist,fn1)

                # Calculate new f_n given new ht4 value
                fn = (ht4 - ht3)/(etab*hpr - ht4)

                # Stopping infinite loop
                if n == nmax:
                    print('WARNING- Max iterations met')
            return fn,ht4
        else:
            ht3 = cubic_spline_interpolation(TC,hf0,Tt3) 
            f0 = flist[0]

            # Initial guess of f
            fn1 = flist[1] 

            # Interpolated values of h at Tt4 for each value of f
            hlist = [cubic_spline_interpolation(TC,hf0,Tt4),cubic_spline_interpolation(TC,hf169,Tt4),cubic_spline_interpolation(TC,hf338,Tt4),cubic_spline_interpolation(TC,hf507,Tt4),cubic_spline_interpolation(TC,hf676,Tt4)]
            
            # Initial guess of ht4 given initial guess of f
            ht4 = hlist[1]
            fn = (ht4-ht3)/(etab*hpr - ht4)
            n = 0

            # Iterating until the stopping conditions are met
            while np.abs(fn-fn1)>tol and n<nmax:
                if pr:
                    print(f"f = {fn}, ht4 = {ht4}, n = {n}")

                # Set new f_{n-1}
                fn1 = fn
                # Increase n to prevent infinite looping
                n+=1

                # Calculate ht4 at new value of f_{n-1}
                ht4 = cubic_spline_interpolation(flist,hlist,fn1)

                # Calculate new f_n given new ht4 value
                fn = (ht4 - ht3)/(etab*hpr - ht4)

                # Stopping infinite loop
                if n == nmax:
                    print('WARNING- Max iterations met')
            return fn,ht4
    