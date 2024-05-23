def Notch_Fliege(f0, Q = 10, margin = 2, RQ = None, R0 = None, C0 = None, cap = None, res = None, visible='Off'):
    ### Required libraries are imported
    import si_prefix
    import numpy as np
    from tabulate import tabulate
    
    #
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    
    ### Initialize capacitor values
    if cap is None:
        # Default capacitor values
        cap_p = np.array([0.5, 0.75, 1, 1.2, 1.3, 1.5, 1.6, 1.8, 2, 2.2, 2.4, 2.5, 2.7, 3, 3.3, 3.6, 3.9, 4.3, 4.7, \
         5, 6, 6.2, 6.8, 7, 7.5, 8, 8.2, 9, 10, 11, 12, 15, 16, 18, 20, 22, 24, 27, 30, 33, 36, 39,\
         43, 47, 51, 56, 62, 68, 82, 100, 110, 120, 150, 180, 200, 220, 330, 470, 680, 820])*(10**-12)

        cap_n = np.array([1, 1.5, 2.2, 3.3, 4.7, 5.6, 6.8, 8.2, 10, 15, 22, 33, 47, 68, 82, 100, 220, 470, 680, 820])\
                         *(10**-9)

        cap_u = np.array([1, 1.5, 2.2, 3.3, 4.7, 5.6])*(10**-6)

        cap_library = np.concatenate((cap_p,cap_n,cap_u),axis=0)
        
    elif type(cap) is np.ndarray:
        # Assign custom capacitor values
        cap_library = cap
        
    else:
        raise TypeError("Capacitor values must be np.ndarray data type!")
    

    ### Initialize resistor values
    if res is None:
        # Default resistor values
        res_0 = np.array([1, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2, 2.2, 2.4, 2.7, 3, 3.3, 3.6,\
                          3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1])
        res_1 = res_0*10
        res_2 = res_1*10
        res_3 = res_2*10
        res_4 = res_3*10
        res_5 = res_4*10
        res_6 = res_5*10
        res_library = np.concatenate(([0],res_0,res_1,res_2,res_3,res_4,res_5,res_6,[10**7]),axis=0)
        
    elif type(res) is np.ndarray:
        # Assign custom resistor values
        res_library = res
    
    else:
        raise TypeError("Resistor values must be np.ndarray data type!")
    
    
    ### Q value 
    Q_ideal = Q
    Q_margin = (margin/100)*Q_ideal
    
    ### Assign cut-off frequency
    f0_ideal = f0
    
    ### Start R and C value estimation
    R0C0 = 2*np.pi*f0_ideal
    C0_value = []
    R0_value, RQ_value = [], []
    Q_value, f0_value = [], []
    for c_temp in cap_library:
        # Start assigning C0 capacitor value
        C0_candid = c_temp
        
        # Find the closest C0 values from capacitor library
        simple_C0 = si_prefix.split(C0_candid)
        if (simple_C0[1] == -6 and simple_C0[0]>= 6) or simple_C0[1] > -6:
            continue
        else:
            C0_candid = find_nearest(cap_library,C0_candid)
            # print("C0 sectim!")
        # Estimate ideal R0
        R0_ideal = (1/R0C0) / C0_candid
        # RQ_ideal = 2 * R0_ideal * Q_ideal
        # Find the closest R0 values from resistor library
        simple_R0 = si_prefix.split(R0_ideal)
        if (simple_R0[1] == 6 and simple_R0[0] >= 10) or simple_R0[1] > 6:
            continue
        else:
            R0_candid_base = np.around(np.array(simple_R0[0]),decimals=1)
            R0_candid_power = (10**simple_R0[1])
            R0_candid = find_nearest(res_library,R0_candid_base*R0_candid_power)
            # print("R0 sectim!")

            RQ_ideal = 2 * R0_candid * Q_ideal
            simple_RQ = si_prefix.split(RQ_ideal)
            if (simple_RQ[1] == 6 and simple_RQ[0] >= 10) or simple_RQ[1] > 6:
                continue
            else:
                RQ_candid_base = np.around(np.array(simple_RQ[0]),decimals=1)
                RQ_candid_power = (10**simple_RQ[1])        
                RQ_candid = find_nearest(res_library,RQ_candid_base*RQ_candid_power)

            f0_candid = 1/(2 * np.pi * R0_candid * C0_candid)
            Q_candid = RQ_candid/(2 * R0_candid)

            # Check if actual Q values is between decided margin values
            if np.abs(Q_candid-Q_ideal) <= Q_margin:
                # Append capacitor and resistor values to empty array
                C0_value.append(si_prefix.si_format(C0_candid))
                R0_value.append(si_prefix.si_format(R0_candid))
                RQ_value.append(si_prefix.si_format(RQ_candid))
                
                Q_value.append(Q_candid)
                f0_value.append(f0_candid)
    
    ### Create dictionary for estimated values 
    selected_values = {'C0':C0_value, 'R0':R0_value,'RQ':RQ_value,\
        'Q factor':Q_value,'Frequency Center':f0_value}
    
    ### Return all possible R/C combinations
    if visible.casefold() == 'on':
        print(tabulate(selected_values,headers='keys'))
        return selected_values
    else:
        return selected_values
