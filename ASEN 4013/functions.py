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
                    print(f"    f = {fn}, ht4 = {ht4}, n = {n}")

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
            if pr:
                print(f"    f = {fn}, ht4 = {ht4}, n = {n}")
            return fn,ht3,ht4
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
                    print(f"    f = {fn}, ht4 = {ht4}, n = {n}")

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
            if pr:
                    print(f"    f = {fn}, ht4 = {ht4}, n = {n}")
            return fn,ht3,ht4
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
            return fn,ht3,ht4
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
            return fn,ht3,ht4
        


def pca(givens:dict,model:str = 'TURBOJET',SI:bool =True)->dict:
    """  
    This function performs Parametric Cycle Analysis for a real engine using the algorithms 
    described in Chapter 7 of Mattingly. 
    
    Args:
        givens: A dict containing all the necessary values, given or assumed. See Mattingly for required info
        model: String dictating which algorithm to follow (TURBOJET or TURBOFAN)
        SI: bool for what type of units, defaults to True (meaning SI). Also assumes then in kJ/kg per table entries

    Returns:
        solution: dict containing the helpful values to return
    """
    #for turbofan dev
    #M0,T0,gac,cpc,gat,cpt,hpr,pid_max,pib,pin,ec,et,etab,etam,P0P9,Tt4,pic
#    givens_fan = {"M0":0.8,"T0":390,"gac":1.4,"cpc":0.24,"gat":1.33,"cpt":0.276,"hpr":18400,"pid_max":0.99,
#           "pib":0.96,"pin":0.99,"ec":0.9,"et":0.9,"etab":0.99,"etam":0.99,"P0P9":0.5,"Tt4":2420,"pic":10,"pifn":0.99,
#           "ef":0.89,"pif":1.65,"alpha":7}

    if model.upper() == 'TURBOJET':
        M0,T0,gac,cpc,gat,cpt,hpr,pid_max,pib,pin,ec,et,etab,etam,P0P9,Tt4,pic = givens.items()
        M0 = M0[1]
        T0 = T0[1]
        gac = gac[1]
        cpc = cpc[1]
        gat = gat[1]
        cpt = cpt[1]
        hpr = hpr[1]
        pid_max = pid_max[1]
        pib = pib[1]
        pin = pin[1]
        ec = ec[1]
        et = et[1]
        etab = etab[1]
        etam = etam[1]
        P0P9 = P0P9[1]
        Tt4 = Tt4[1]
        pic = pic[1]

        # Defining gc
        if SI:
            gc = 1
        else:
            gc = 32.174

        print('Freestream and Ram Properties')
        Rc = (gac - 1)/gac*cpc
        Rt = (gat - 1)/gat*cpt
        if SI:
            a0 = np.sqrt(gac*Rc*gc*T0*1000)
        V0 = a0*M0
        taur = 1 + (gac-1)/2*M0**2
        pir = taur**(gac/(gac-1))
        # Finding eta_r
        if M0<1:
            etar = 1
        elif M0>1 and M0<5:
            etar = 1 - 0.075*(M0 - 1)**1.35
        else:
            etar = 800/(M0**4 + 935)
        print(f'    R_c = {round(Rc,4)} kJ/kgK, R_t = {round(Rt,4)} kJ/kgK')
        print(f'    a_0 = {round(a0,4)} m/s, V_0 = {round(V0,4)} m/s')
        print(f'    tau_r = {round(taur,4)}, pi_r = {round(pir,4)}, eta_r = {round(etar,4)}')

        print('Inlet Diffuser')
        pid = pid_max*etar
        taul = cpt*Tt4/(cpc*T0)
        print(f'    pi_d = {round(pid,4)}, tau_lambda = {round(taul,4)}')

        print('Compressor')
        tauc = pic**((gac - 1)/(gac*ec))
        etac = (pic**((gac-1)/gac)-1)/(tauc - 1)
        Tt3 = T0*taur*tauc
        print(f'    tau_c = {round(tauc,4)}, eta_c = {round(etac,4)}, Tt3 = {round(Tt3,4)}')

        print('Iterating to find ht4 and f')
        f,ht3,ht4 = find_f(Tt3,Tt4,etab,hpr,pr = True)
        print(f'    Final Values: f = {round(float(f),6)}, ht4 = {round(float(ht4),4)} kJ/kg')

        print('Turbine Properties')
        taut = 1 - 1/(etam*(1+f))*taur/taul*(tauc-1)
        pit = taut**(gat/((gat-1)*et))
        etat = (1 - taut)/(1-taut**(1/et))
        print(f'    tau_t = {round(taut,4)}, pi_t = {round(pit,4)}, eta_t = {round(etat, 4)}')

        print('Exit Properties')
        Pt9P9 = P0P9*pir*pid*pic*pib*pit*pin
        M9 = np.sqrt(2/(gat-1)*(Pt9P9**((gat-1)/gat)-1))
        T9T0 = Tt4*taut/T0/(Pt9P9**((gat-1)/gat))
        T9 = T9T0*T0
        V9a0 = M9*np.sqrt(gat*Rt*T9/(gac*Rc*T0))
        print(f'    Pt9/P9 = {round(Pt9P9,4)}, T9 = {round(T9,4)} K, T9/T0 = {round(T9T0,4)}')
        print(f'    M9 = {round(M9,4)}, V9/a0 = {round(V9a0,4)}')

        print('Performance Properties')
        Fm0 = a0/gc*((1+f)*V9a0 - M0 + (1+f)*Rt*T9T0/(Rc*V9a0)*(1-P0P9)/gac)
        S = f/Fm0*1e6
        etaTH = a0**2*((1+f)*V9a0**2 - M0**2)/(2*gc*f*hpr)/1000
        etap = 2*gc*V0*Fm0/(a0**2*((1+f)*V9a0**2 - M0**2))
        etaO = etaTH*etap
        print(f'    F/m0 = {round(Fm0,4)} N/kg/s, S = {round(S,4)} mg/s/kN')
        print(f'    eta_TH = {round(etaTH,4)}, eta_p = {round(etap,4)}, eta_O = {round(etaO,4)}')

        solution = {"Rc":Rc,"Rt":Rt,"a0":float(a0),"V0":float(V0),"taur":taur,"pir":pir,"etar":etar,"pid":pid,
                    "taul":taul,"tauc":tauc,"etac":etac,"Tt3":Tt3,"f":float(f),"ht3":float(ht3),"taut":float(taut),
                    "pit":float(pit),"etat":float(etat),"Pt9P9":float(Pt9P9),"M9":float(M9),"T9T0":float(T9T0),"V9a0":float(V9a0),"Fm0":float(Fm0),
                    "S":float(S),"etaTH":float(etaTH),"etap":float(etap),"etaO":float(etaO)}
        return solution
    elif model.upper() == 'TURBOFAN':
        print('WARNING: TURBOFAN NOT COMPLETED')
        M0,T0,gac,cpc,gat,cpt,hpr,pid_max,pib,pin,ec,et,etab,etam,P0P9,Tt4,pic,pifn,ef,pif,alpha = givens.items()
        M0 = M0[1]
        T0 = T0[1]
        gac = gac[1]
        cpc = cpc[1]
        gat = gat[1]
        cpt = cpt[1]
        hpr = hpr[1]
        pid_max = pid_max[1]
        pib = pib[1]
        pin = pin[1]
        ec = ec[1]
        et = et[1]
        etab = etab[1]
        etam = etam[1]
        P0P9 = P0P9[1]
        Tt4 = Tt4[1]
        pic = pic[1]

        # Defining gc
        if SI:
            gc = 1
        else:
            gc = 32.174

        print('Freestream and Ram Properties')
        Rc = (gac - 1)/gac*cpc
        Rt = (gat - 1)/gat*cpt
        if SI:
            a0 = np.sqrt(gac*Rc*gc*T0*1000)
        V0 = a0*M0
        taur = 1 + (gac-1)/2*M0**2
        pir = taur**(gac/(gac-1))
        # Finding eta_r
        if M0<1:
            etar = 1
        elif M0>1 and M0<5:
            etar = 1 - 0.075*(M0 - 1)**1.35
        else:
            etar = 800/(M0**4 + 935)
        print(f'    R_c = {round(Rc,4)} kJ/kgK, R_t = {round(Rt,4)} kJ/kgK')
        print(f'    a_0 = {round(a0,4)} m/s, V_0 = {round(V0,4)} m/s')
        print(f'    tau_r = {round(taur,4)}, pi_r = {round(pir,4)}, eta_r = {round(etar,4)}')

        print('Inlet Diffuser')
        pid = pid_max*etar
        taul = cpt*Tt4/(cpc*T0)
        print(f'    pi_d = {round(pid,4)}, tau_lambda = {round(taul,4)}')

        print('Compressor')
        tauc = pic**((gac - 1)/(gac*ec))
        etac = (pic**((gac-1)/gac)-1)/(tauc - 1)
        Tt3 = T0*taur*tauc
        print(f'    tau_c = {round(tauc,4)}, eta_c = {round(etac,4)}, Tt3 = {round(Tt3,4)}')

        print('Iterating to find ht4 and f')
        f,ht3,ht4 = find_f(Tt3,Tt4,etab,hpr,pr = True)
        print(f'    Final Values: f = {round(float(f),6)}, ht4 = {round(float(ht4),4)} kJ/kg')

        print('Turbine Properties')
        taut = 1 - 1/(etam*(1+f))*taur/taul*(tauc-1)
        pit = taut**(gat/((gat-1)*et))
        etat = (1 - taut)/(1-taut**(1/et))
        print(f'    tau_t = {round(taut,4)}, pi_t = {round(pit,4)}, eta_t = {round(etat, 4)}')

        print('Exit Properties')
        Pt9P9 = P0P9*pir*pid*pic*pib*pit*pin
        M9 = np.sqrt(2/(gat-1)*(Pt9P9**((gat-1)/gat)-1))
        T9T0 = Tt4*taut/T0/(Pt9P9**((gat-1)/gat))
        T9 = T9T0*T0
        V9a0 = M9*np.sqrt(gat*Rt*T9/(gac*Rc*T0))
        print(f'    Pt9/P9 = {round(Pt9P9,4)}, T9 = {round(T9,4)} K, T9/T0 = {round(T9T0,4)}')
        print(f'    M9 = {round(M9,4)}, V9/a0 = {round(V9a0,4)}')

        print('Performance Properties')
        Fm0 = a0/gc*((1+f)*V9a0 - M0 + (1+f)*Rt*T9T0/(Rc*V9a0)*(1-P0P9)/gac)
        S = f/Fm0*1e6
        etaTH = a0**2*((1+f)*V9a0**2 - M0**2)/(2*gc*f*hpr)/1000
        etap = 2*gc*V0*Fm0/(a0**2*((1+f)*V9a0**2 - M0**2))
        etaO = etaTH*etap
        print(f'    F/m0 = {round(Fm0,4)} N/kg/s, S = {round(S,4)} mg/s/kN')
        print(f'    eta_TH = {round(etaTH,4)}, eta_p = {round(etap,4)}, eta_O = {round(etaO,4)}')

        solution = {"Rc":Rc,"Rt":Rt,"a0":float(a0),"V0":float(V0),"taur":taur,"pir":pir,"etar":etar,"pid":pid,
                    "taul":taul,"tauc":tauc,"etac":etac,"Tt3":Tt3,"f":float(f),"ht3":float(ht3),"taut":float(taut),
                    "pit":float(pit),"etat":float(etat),"Pt9P9":float(Pt9P9),"M9":float(M9),"T9T0":float(T9T0),"V9a0":float(V9a0),"Fm0":float(Fm0),
                    "S":float(S),"etaTH":float(etaTH),"etap":float(etap),"etaO":float(etaO)}
        return solution
    