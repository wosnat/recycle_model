import numpy as np
    
def calc_csat(T = None,S = None,pco2eq = None,pt = None,sit = None,ta = None): 
    ##########################################################
    """
        # CALCULATE Csat, EQUILIBRIUM DIC FOR GIVEN pCO2, T, Alk etc...
    # Efficient solver following Follows et al (2005)
    #       pco2eq = atmospheric reference pCO2 level (atmospheres)
    #                for which to find equilibrium dic, csat
    #       csat = equilibrium total inorganic carbon (mol/m^3)
    #             where 1 T = 1 metric ton = 1000 kg
    #       ta  = total alkalinity (eq/m^3)
    #       pt  = inorganic phosphate (mol/^3)
    #       sit = inorganic silicate (mol/^3)
    #       T   = temperature (degrees K)
    #       S   = salinity (PSU)
    #		hg  = first guess of [H+] (10e-8 is a good cold start value)
    #       convert mol kg-1, mol m-3
    #       permil = 1.0 / 1024.5
    """    
    # convert from  mol/m^3 to mol/kg - nominal density of water set at 1024.5
    permil = 1.0 / 1024.5
    pt = pt * permil
    #   dic = dic*permil;
    sit = sit * permil
    ta = ta * permil
    # set first guess for [H]
    hg = 1e-08
    # estimated concentration of borate based on salinity
    scl = S / 1.80655
    bt = 0.000232 * scl / 10.811
    # some definitions ...
    S2 = S * S
    sqrtS = S ** 0.5
    invT = 1.0 / T
    T1 = T / 100.0
    # Coefficient algorithms as used in OCMIP2 protocols
# K1, K2 Millero (1995) using Mehrbach data
    k1 = 10 ** (- 1 * (3670.7 * invT - 62.008 + 9.7944 * np.log(T) - 0.0118 * S + 0.000116 * S2))
    k2 = 10 ** (- 1 * (1394.7 * invT + 4.777 - 0.0184 * S + 0.000118 * S2))
    # K1p, K2p, K3p  DOE (1994)
    k1p = np.exp(- 4576.752 * invT + 115.525 - 18.453 * np.log(T) + (- 106.736 * invT + 0.69171) * sqrtS + (- 0.65643 * invT - 0.01844) * S)
    k2p = np.exp(- 8814.715 * invT + 172.0883 - 27.927 * np.log(T) + (- 160.34 * invT + 1.3566) * sqrtS + (0.37335 * invT - 0.05778) * S)
    k3p = np.exp(- 3070.75 * invT - 18.141 + (17.27039 * invT + 2.81197) * sqrtS + (- 44.99486 * invT - 0.09984) * S)
    # Kb, Millero (1995) using data from Dickson
    kb = np.exp((- 8966.9 - 2890.53 * sqrtS - 77.942 * S + 1.728 * S ** 1.5 - 0.0996 * S2) * invT + (148.0248 + 137.1942 * sqrtS + 1.62142 * S) + (- 24.4344 - 25.085 * sqrtS - 0.2474 * S) * np.log(T) + 0.053105 * T * sqrtS)
    # Kw, Millero (1995)
    kw = np.exp(- 13847.26 * invT + 148.9652 - 23.6521 * np.log(T) + (118.67 * invT - 5.977 + 1.0495 * np.log(T)) * sqrtS - 0.01615 * S)
    # Ksi, Millero (1995)
    I = (19.924 * S) / (1000 - 1.005 * S)
    ksi = np.exp(- 8904.2 * invT + 117.385 - 19.334 * np.log(T) + (- 458.79 * invT + 3.5913) * (I ** 0.5) + (188.74 * invT - 1.5998) * I + (- 12.1652 * invT + 0.07871) * (I * I) + np.log(1.0 - 0.001005 * S))
    # fugacity, Weiss and Price, Marine Chem, 8, 347 (1990)
    ff = np.exp(- 162.8301 + 218.2968 / (T1) + 90.9241 * np.log(T1) - 1.47696 * (T1 * T1) + S * (0.025695 - 0.025225 * (T1) + 0.0049867 * (T1 * T1)))
    # First guess of [H+]: from last timestep *OR* fixed for cold start
# --- here iterate for accurate solution
    for i in range(10):
        # estimate contributions to total alk from borate, silicate, phosphate
        bohg = (bt * kb) / (hg + kb)
        siooh3g = (sit * ksi) / (ksi + hg)
        denom = (hg * hg * hg) + (k1p * hg * hg) + (k1p * k2p * hg) + (k1p * k2p * k3p)
        h3po4g = (pt * hg * hg * hg) / denom
        h2po4g = (pt * k1p * hg * hg) / denom
        hpo4g = (pt * k1p * k2p * hg) / denom
        po4g = (pt * k1p * k2p * k3p) / denom
        # estimate carbonate alkalinity
        cag = ta - bohg - (kw / hg) + hg - hpo4g - 2 * po4g + h3po4g - siooh3g
        # estimate hydrogen ion conc
        stuff = k1 * ff * pco2eq / cag
        stuff2 = stuff * (stuff + 8.0 * k2)
        stuff3 = stuff2 ** 0.5
        Hnew = 0.5 * (stuff + stuff3)
        hg = Hnew
        #       pH = -log10(Hnew);
#       fprintf('pH #14.10e\n',pH);
    
    # evaluate csat, equilibrium DIC concentration
    csat = pco2eq * ff * (1.0 + (k1 / Hnew) + (k1 * k2 / (Hnew * Hnew)))
    # Output pt, sit, ta and csat in mol/m^3
    pt = pt / permil
    sit = sit / permil
    ta = ta / permil
    csat = csat / permil
    # calc final pH
    pH = - np.log10(Hnew)
    # some diagnostics.........
#fprintf('T #10.4f\n',T);
#fprintf('S #10.4f\n',S);
#fprintf('Alk #10.4f\n',ta);
#fprintf('bt #10.4e\n',bt);
#fprintf('k1 #10.4e\n',k1);
#fprintf('k2 #10.4e\n',k2);
#fprintf('k1p #10.4e\n',k1p);
#fprintf('k2p #10.4e\n',k2p);
#fprintf('k3p #10.4e\n',k3p);
#fprintf('kw #10.4e\n',kw);
#fprintf('kb #10.4e\n',kb);
#fprintf('pH #10.4f\n',pH);
#fprintf('csat #10.4f\n',csat);
    
    return csat
    
