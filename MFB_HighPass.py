def MFB_HighPass(cutoff, gain = 1, Q = 0.707, margin = 2, filter_type = None, cap = None, res = None, visible='Off'):
    ### Required libraries are imported
    import si_prefix
    import numpy as np
    from tabulate import tabulate
    
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
    
    
    ### Filter type and Q value 
    # Available filter types are 'bessel' , 'linear' , 'butter'
    if filter_type is None:
        # No filter type is specified
        # Look for Q value
        Q_ideal = Q
        Q_margin = (margin/100)*Q_ideal
    elif type(filter_type) is not str:
        raise TypeError("filter type must be string!")
    else:
        try:
            filters = {'bessel': 0.577, 'linear': 0.644, 'butter': 0.707}
            Q_ideal = filters[filter_type.casefold()]
            Q_margin = (margin/100)*Q_ideal
        except KeyError:
            print('Filter type must be "Bessel" , "Linear" or "Butter" ')
            raise
    
    ### Assign ideal gain value
    gain_ideal = gain
    
    ### Assign cut-off frequency
    fc_ideal = cutoff
    
    ### Start R and C value estimation
    wc = 2*np.pi*fc_ideal
    C1_value, C2_value, C3_value = [], [], []
    R1_value, R2_value = [], []
    Q_value, Gain_value, fc_value = [], [], []
    for c_temp in cap_library:
        # Start assigning C3 capacitor value
        C3_candid = c_temp

        # Find ideal C1 and C2 values
        C1_ideal = c_temp*gain_ideal
        C2_ideal = c_temp*gain_ideal
        
        # Find the closest C1 and C2 values from capacitor library
        simple_C1 = si_prefix.split(C1_ideal)
        simple_C2 = si_prefix.split(C2_ideal)
        if (simple_C1[1] == -6 and simple_C1[0]>= 6) or simple_C1[1] > -6:
            continue
        else:
            C1_candid = find_nearest(cap_library,C1_ideal)
            C2_candid = find_nearest(cap_library,C2_ideal)
        
        # Estimate ideal R1 and R2 values
        gain_actual = C1_candid/C3_candid
        R1_ideal = (1/Q_ideal)/(wc*C1_candid*(2+1/gain_actual))
        R2_ideal = (gain_actual*(2+1/gain_actual))/(wc*C1_candid*(1/Q_ideal))
        
        # Find the closest R1 and R2 values from resistor library
        simple_R1 = si_prefix.split(R1_ideal)
        simple_R2 = si_prefix.split(R2_ideal)
        if (simple_R1[1] == 6 and simple_R1[0] >= 10) or simple_R1[1] > 6:
            continue
        elif (simple_R2[1] == 6 and simple_R2[0] >= 10) or simple_R2[1] > 6:
            continue
        else:
            R1_candid_base = np.around(np.array(simple_R1[0]),decimals=1)
            R1_candid_power = (10**simple_R1[1])
            R1_candid = find_nearest(res_library,R1_candid_base*R1_candid_power)

            R2_candid_base = np.around(np.array(simple_R2[0]),decimals=1)
            R2_candid_power = (10**simple_R2[1])        
            R2_candid = find_nearest(res_library,R2_candid_base*R2_candid_power)

            w0 = 1/(R1_candid*R2_candid*C2_candid*C3_candid)

            Q_candid = np.sqrt(R2_candid*C2_candid*C3_candid)/\
                ((C1_candid+C2_candid+C3_candid)*np.sqrt(R1_candid))
            fc_candid = np.sqrt(1/(R1_candid*R2_candid*C2_candid*C3_candid))/(2*np.pi)
        
        # Check if actual Q values is between decided margin values
        if np.abs(Q_candid-Q_ideal) <= Q_margin:
            # Append capacitor and resistor values to empty array
            C1_value.append(si_prefix.si_format(C1_candid))
            C2_value.append(si_prefix.si_format(C2_candid))
            C3_value.append(si_prefix.si_format(C3_candid))
            R1_value.append(si_prefix.si_format(R1_candid))
            R2_value.append(si_prefix.si_format(R2_candid))
            
            Q_value.append(Q_candid)
            Gain_value.append(gain_actual)
            fc_value.append(fc_candid)
    
    ### Create dictionary for estimated values 
    selected_values = {'C1':C1_value, 'C2':C2_value,'C3':C3_value,'R1':R1_value,'R2':R2_value,\
        'Q factor':Q_value,'Cut-off':fc_value,'Gain':Gain_value}
    
    ### Return all possible R/C combinations
    if visible.casefold() == 'on':
        print(tabulate(selected_values,headers='keys'))
        return selected_values
    else:
        return selected_values
