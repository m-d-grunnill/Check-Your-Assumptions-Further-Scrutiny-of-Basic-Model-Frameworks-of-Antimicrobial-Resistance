# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:09:58 2019
Updated on Tue May 07 2019
Update on Wed Nov 11 2021

@author: Martin.Grunnill

The first cell of this module has four functions that setup Single-Strain, the Superinfection, the 
Exclusive infection, the replacement infection and the bi-directional conversion 
models (Spicknall et al 2013) in PyGOM.

The second cell of this module has three functions that determine the endemic equilibria for the 
Exclusive infection, the replacement infection and the bi-directional conversion 
models as outlined in Spicknall et al (2013). 
"""

#%% PyGOM Models

import pygom

def AMR_single_strain_intia():
    """This function sets up the single strain model as Outlined in 
        Spicknall et al (2013)."""
    states = ['S_W','I_W','S_Z','I_Z']
    params = ['N','beta_W','beta_Z','gamma','gamma_T','epsilon']
    
    S_W_to_I_W = pygom.Transition (origin='S_W',
                                   equation ='(beta_W/N)*I_W*S_W',
                                   transition_type = 'T',
                                   destination ='I_W')
    
    S_Z_to_I_Z = pygom.Transition (origin='S_Z',
                                   equation ='(beta_Z/N)*I_Z*S_Z',
                                   transition_type = 'T',
                                   destination ='I_Z')

    I_W_to_S_W = pygom.Transition (origin='I_W',
                                   equation ='I_W*epsilon*gamma_T + I_W*gamma',
                                   transition_type = 'T',
                                   destination ='S_W')
    
    I_Z_to_S_Z = pygom.Transition (origin='I_Z',
                                  equation ='gamma*I_Z',
                                  transition_type = 'T',
                                  destination ='S_Z')

    model = pygom.SimulateOde (states ,
                               params ,
                               transition = [S_W_to_I_W,
                                             S_Z_to_I_Z,
                                             I_W_to_S_W,                                             
                                             I_Z_to_S_Z])
    return(model)
    

def AMR_superinf_intia():
    """This function setsup Superinfection model as Outlined in 
        Spicknall et al (2013). This includes Single-Strain AMR model , as the 
        difference between these models is the setting of a parameter to 0 or 1."""
    
    states = ['S','I_W','I_Z','I_WZ']
    params = ['N','beta_W','beta_Z','gamma','gamma_T','epsilon','q']
    
    S_to_I_W = pygom.Transition (origin='S',
                             equation ='(beta_W/N)*(I_W+q*I_WZ)*S',
                             transition_type = 'T',
                             destination ='I_W')

    S_to_I_Z = pygom.Transition (origin='S',
                                 equation ='(beta_Z/N)*(I_Z+q*I_WZ)*S',
                                 transition_type = 'T',
                                 destination ='I_Z')
    
    I_W_to_S = pygom.Transition (origin='I_W',
                                 equation ='gamma*(1-epsilon)*I_W + gamma_T*epsilon*I_W',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_W_to_I_WZ = pygom.Transition (origin='I_W',
                                    equation ='(beta_Z/N)*I_W*(I_Z+q*I_WZ)',
                                    transition_type = 'T',
                                    destination ='I_WZ')
    
    I_Z_to_S = pygom.Transition (origin='I_Z',
                                 equation ='gamma*I_Z',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_Z_to_I_WZ = pygom.Transition (origin='I_Z',
                                    equation ='(beta_W/N)*I_Z*(I_W+q*I_WZ)',
                                    transition_type = 'T',
                                    destination ='I_WZ')
    
    I_WZ_to_I_W = pygom.Transition (origin='I_WZ',
                                    equation ='gamma*I_WZ',
                                    transition_type = 'T',
                                    destination ='I_W')
    
    I_WZ_to_I_Z = pygom.Transition (origin='I_WZ',
                                    equation ='(1-epsilon)*gamma*I_WZ + epsilon*gamma_T*I_WZ',
                                    transition_type = 'T',
                                    destination ='I_Z')
    
    model = pygom.SimulateOde (states ,
                               params ,
                               transition = [S_to_I_W,
                                             S_to_I_Z,
                                             I_W_to_S,
                                             I_W_to_I_WZ,
                                             I_Z_to_S,
                                             I_Z_to_I_WZ,
                                             I_WZ_to_I_W,
                                             I_WZ_to_I_Z])
    return(model)


def AMR_full_co_inf_intia():
    """This function sets up Bi & Unidirectional conversions AMR
    models as Outlined in  Spicknall et al (2013). Conversion from sensitive to
    resistant infections occurs when rho is >0. Conversion from resistant to
    sensitive infections occurs when phi is >0. """
    # Create coinfection model
    states = ['S', 'I_W', 'I_Z', 'I_WZ', 'I_WW', 'I_ZZ']
    params = ['N', 'beta_W', 'beta_Z', 'gamma', 'gamma_T', 'epsilon', 'q']

    foi_W = '(beta_W/N)*(I_W+q*(I_WZ+I_WW+I_WW))'
    foi_Z = '(beta_Z/N)*(I_Z+q*(I_WZ+I_ZZ+I_ZZ))'

    S_to_I_W = pygom.Transition(origin='S',
                                equation='S*' + foi_W,
                                transition_type='T',
                                destination='I_W')

    S_to_I_Z = pygom.Transition(origin='S',
                                equation='S*' + foi_Z,
                                transition_type='T',
                                destination='I_Z')

    I_W_to_S = pygom.Transition(origin='I_W',
                                equation='gamma*(1-epsilon)*I_W + gamma_T*epsilon*I_W',
                                transition_type='T',
                                destination='S')

    I_Z_to_S = pygom.Transition(origin='I_Z',
                                equation='gamma*I_Z',
                                transition_type='T',
                                destination='S')

    I_W_to_I_WZ = pygom.Transition(origin='I_W',
                                   equation='I_W*' + foi_Z,
                                   transition_type='T',
                                   destination='I_WZ')

    I_Z_to_I_WZ = pygom.Transition(origin='I_Z',
                                   equation='I_Z*' + foi_W,
                                   transition_type='T',
                                   destination='I_WZ')

    I_WZ_to_I_W = pygom.Transition(origin='I_WZ',
                                   equation='gamma*I_WZ',
                                   transition_type='T',
                                   destination='I_W')

    I_WZ_to_I_Z = pygom.Transition(origin='I_WZ',
                                   equation='((1-epsilon)*gamma + epsilon*gamma_T)*I_WZ',
                                   transition_type='T',
                                   destination='I_Z')

    I_W_to_I_WW = pygom.Transition(origin='I_W',
                                   equation='I_W*' + foi_W,
                                   transition_type='T',
                                   destination='I_WW')

    I_WW_to_S = pygom.Transition(origin='I_WW',
                                 equation='((1-epsilon)*gamma + gamma_T*epsilon)*I_WW',
                                 transition_type='T',
                                 destination='S')

    I_Z_to_I_ZZ = pygom.Transition(origin='I_Z',
                                   equation='I_Z*' + foi_Z,
                                   transition_type='T',
                                   destination='I_ZZ')

    I_ZZ_to_S = pygom.Transition(origin='I_ZZ',
                                 equation='I_ZZ*gamma',
                                 transition_type='T',
                                 destination='S')

    model = pygom.SimulateOde(states,
                              params,
                              transition=[S_to_I_W,
                                          S_to_I_Z,
                                          I_W_to_S,
                                          I_W_to_I_WZ,
                                          I_Z_to_S,
                                          I_Z_to_I_WZ,
                                          I_WZ_to_I_W,
                                          I_WZ_to_I_Z,
                                          I_W_to_I_WW,
                                          I_WW_to_S,
                                          I_Z_to_I_ZZ,
                                          I_ZZ_to_S])

    return model
    

def AMR_exclus_inf_intia():
    """This function setsup Exclusive infection model as Outlined in 
        Spicknall et al (2013)."""
     
    states = ['S','I_W','I_Z']
    params = ['N','beta_W','beta_Z','gamma','gamma_T','epsilon']
    
    S_to_I_W = pygom.Transition (origin='S',
                             equation ='I_W*S*(beta_W/N) + I_W*epsilon*gamma',
                             transition_type = 'T',
                             destination ='I_W')

    S_to_I_Z = pygom.Transition (origin='S',
                                 equation ='I_Z*S*(beta_Z/N)',
                                 transition_type = 'T',
                                 destination ='I_Z')
    
    I_W_to_S = pygom.Transition (origin='I_W',
                                 equation ='I_W*epsilon*gamma_T + I_W*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_Z_to_S = pygom.Transition (origin='I_Z',
                                 equation ='I_Z*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    model = pygom.SimulateOde (states ,
                               params ,
                               transition = [S_to_I_W,
                                             S_to_I_Z,
                                             I_W_to_S,
                                             I_Z_to_S])
    
    return(model)

def AMR_replace_inf_intia():
    """This function setsup replacement infection model as Outlined in 
        Spicknall et al (2013). When ((beta_W/N)-(beta_Z/N)) = 0 
        there is no replacement infection model becomes the exclusive 
        infection model."""
     
    states = ['S','I_W','I_Z']
    params = ['N','beta_W','beta_Z','gamma','gamma_T','epsilon']
    
    S_to_I_W = pygom.Transition (origin='S',
                             equation ='I_W*S*(beta_W/N) + I_W*epsilon*gamma',
                             transition_type = 'T',
                             destination ='I_W')

    S_to_I_Z = pygom.Transition (origin='S',
                                 equation ='I_Z*S*(beta_Z/N)',
                                 transition_type = 'T',
                                 destination ='I_Z')
    
    I_W_to_S = pygom.Transition (origin='I_W',
                                 equation ='I_W*epsilon*gamma_T + I_W*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_W_to_I_Z = pygom.Transition (origin='I_W',
                                   equation ='I_W*I_Z*(beta_Z/N)',
                                   transition_type = 'T',
                                   destination ='I_Z')
    
    I_Z_to_S = pygom.Transition (origin='I_Z',
                                 equation ='I_Z*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_Z_to_I_W = pygom.Transition (origin='I_Z',
                                   equation ='I_W*I_Z*(beta_W/N)',
                                   transition_type = 'T',
                                   destination ='I_W') 
    
    model = pygom.SimulateOde (states ,
                               params ,
                               transition = [S_to_I_W,
                                             S_to_I_Z,
                                             I_W_to_S,
                                             I_W_to_I_Z,
                                             I_Z_to_S,
                                             I_Z_to_I_W])
    
    return(model)


def AMR_bi_and_uni_conversion_inf_intia():
    """This function sets up Bi & Unidirectional conversions AMR 
    models as Outlined in  Spicknall et al (2013). Conversion from sensitive to 
    resistant infections occurs when rho is >0. Conversion from resistant to 
    sensitive infections occurs when phi is >0. """
     
    states = ['S','I_W','I_Z']
    params = ['N','beta_W','beta_Z','gamma','gamma_T','epsilon','rho','phi']
    
   
    S_to_I_W = pygom.Transition (origin='S',
                                 equation ='I_W*S*(beta_W/N) + I_W*epsilon*gamma',
                                 transition_type = 'T',
                                 destination ='I_W')
    
    S_to_I_Z = pygom.Transition (origin='S',
                                 equation ='I_Z*S*(beta_Z/N)',
                                 transition_type = 'T',
                                 destination ='I_Z')
    
    I_W_to_S = pygom.Transition (origin='I_W',
                                 equation ='I_W*epsilon*gamma_T + I_W*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_W_to_I_Z = pygom.Transition (origin='I_W',
                                 equation ='I_W*epsilon*rho + I_Z*epsilon*phi',
                                 transition_type = 'T',
                                 destination ='I_Z')
    
    I_Z_to_S = pygom.Transition (origin='I_Z',
                                 equation ='I_Z*gamma',
                                 transition_type = 'T',
                                 destination ='S')
    
    I_Z_to_I_W = pygom.Transition (origin='I_Z',
                                   equation ='I_Z*phi',
                                   transition_type = 'T',
                                   destination ='I_W')
    
    model = pygom.SimulateOde (states ,
                               params ,
                               transition = [S_to_I_W,
                                             S_to_I_Z,
                                             I_W_to_S,
                                             I_W_to_I_Z,
                                             I_Z_to_S,
                                             I_Z_to_I_W])
    
    return(model)


#%% Function for determining model equilibria stability
import sympy
import math
import copy
from sympy import sqrt

# All parameters and variables need to be setup as 
N, S, I_W, I_Z, beta_W, beta_Z, gamma, gamma_T, epsilon, rho, phi, q =sympy.symbols(
        'N, S I_W I_Z beta_W beta_Z gamma gamma_T epsilon rho phi q')

def round_sf(number, significant):
    '''
    Rounds to a specified number of signicant figures. 
    '''
    return round(number, significant - len(str(number)))

def possible_exclusive_replacement_ee(param_values):
    '''
    Calculates the two two possible endemic equilibria for the Exclusive and Replacement infection models,
    as out lined in Spicknall et al 2013. Reuturning any non-DFE that are
    that are biologically reasonable and locally stable.
    '''
    S_non_DFE_1 = N*(-epsilon*gamma + epsilon * gamma_T + gamma)/beta_W
    S_non_DFE_1 = S_non_DFE_1.subs(param_values)
    EE_1 = {'S': S_non_DFE_1.subs(param_values),
            'I_W': param_values['N']-S_non_DFE_1,
            'I_Z': 0}
    S_non_DFE_2 = N*gamma/beta_Z
    S_non_DFE_2 = S_non_DFE_2.subs(param_values)
    EE_2 = {'S': S_non_DFE_2.subs(param_values),
            'I_W': 0,
            'I_Z': param_values['N']-S_non_DFE_2}
    return {'Endemic Equilibria 1': EE_1, 'Endemic Equilibria 2': EE_2}

def exclus_inf_end_equil(param_values):
    '''
    Calculates the non-disease free equilibria for the exclusive infection model ,
    as out lined in Spicknall et al 2013. Reuturning any non-DFE that are 
    that are biologically reasonable and locally stable.
    '''
    # Note if both strains have R0 values <1 there is no point in using this function.
    
    equil_pops = []
    # Setup Equilibrium populations as an empty list. If any of the endemic 
    # equilibria are found to biologically feasible and stable they are appended to the list.
    
    #Basic Reproductive numbers from Spicknall 2013
    R0_W = beta_W/(gamma*(1-epsilon)+gamma_T*epsilon)
    R0_W = R0_W.subs(param_values)
    R0_Z = beta_Z/gamma
    R0_Z = R0_Z.subs(param_values)
    
    #Formula for the non-Disease Free Equilibria has been worked out in a jupyter notebook
    #and pasted here. NOTE with the first non-DFE the I_W population = N-S (I_Z=0),
    # and with the second non-DFE the I_Z population = N-S (I_W=0).
    # Therefore, only the formule for S are needed (see jupyter notebook).
    S_non_DFE_1 = N*(-epsilon*gamma + epsilon * gamma_T + gamma)/beta_W
    S_non_DFE_1 = S_non_DFE_1.subs(param_values)
    I_W_non_DFE_1 = param_values['N']-S_non_DFE_1
    I_Z_non_DFE_1 = 0
    S_non_DFE_2 = N*gamma/beta_Z
    S_non_DFE_2 = S_non_DFE_2.subs(param_values)
    I_W_non_DFE_2 = 0
    I_Z_non_DFE_2 = param_values['N']-S_non_DFE_2
    
        
    #The Eigen values associated with jacobian matricies for the non-DFE have been 
    #derived in a jupyter notebook and pasted here.
    eigs_J_non_DFE_1 = [-(beta_W*gamma + beta_Z*epsilon*gamma - beta_Z*epsilon*gamma_T - beta_Z*gamma)/beta_W,
                        -I_W*beta_W/N,
                        0]
    eigs_J_non_DFE_2 = [(beta_W*gamma + beta_Z*epsilon*gamma - beta_Z*epsilon*gamma_T - beta_Z*gamma)/beta_Z,
                        -I_Z*beta_Z/N,
                        0]
    
    if math.isnan(S_non_DFE_1) or math.isnan(I_W_non_DFE_1) or math.isnan(I_Z_non_DFE_1) or math.isinf(S_non_DFE_1) or math.isinf(I_W_non_DFE_1) or math.isinf(I_Z_non_DFE_1):
        DFE_1_feasable = False
    else:
        DFE_1_feasable = True
    
    if DFE_1_feasable and round_sf(S_non_DFE_1,2) >= 0 and round_sf(I_W_non_DFE_1,2) >= 0 and round_sf(I_Z_non_DFE_1,2) >= 0:
        DFE_1_feasable = True
    else:
        DFE_1_feasable = False
    
    if math.isnan(S_non_DFE_2) or math.isnan(I_W_non_DFE_2) or math.isnan(I_Z_non_DFE_2) or math.isinf(S_non_DFE_2) or math.isinf(I_W_non_DFE_2) or math.isinf(I_Z_non_DFE_2):
        DFE_2_feasable = False
    else:
        DFE_2_feasable = True
    
    if DFE_2_feasable and round_sf(S_non_DFE_2,2) >= 0 and round_sf(I_W_non_DFE_2,2) >= 0 and round_sf(I_Z_non_DFE_2,2) >= 0:
        DFE_2_feasable = True
    else:
        DFE_2_feasable = False
    
    #Determine if R0_W is >=1, R0_Z<1 and the first DFE is biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)<1 and DFE_1_feasable:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_Z is >=1, R0_W<1 and the second DFE is biologically feasable.
    if R0_Z.subs(param_values)>=1 and R0_W.subs(param_values)<1 and DFE_2_feasable:   
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the first DFE is biologically feasable but the second is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable == False:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the second DFE is biologically feasable but the first is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable  == False and DFE_2_feasable:
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
        
    #Determine if R0_W is >=1, R0_Z>=1 and both DFEs are biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable:   
        # If all the values for non-DFE_1 are >= 0, evaluate eigenavalues 
        # associated with the Jacobian matrix of non-DFE_1.
        # Need to be able to subsititute in equilibria formula, as well as the param values:
        vals_to_subs = copy.deepcopy(param_values)
        vals_to_subs['S'] = S_non_DFE_1
        vals_to_subs['I_W'] = I_W_non_DFE_1
        vals_to_subs['I_Z'] = I_Z_non_DFE_1
        if eigs_J_non_DFE_1[0].subs(vals_to_subs) <=0 and eigs_J_non_DFE_1[1].subs(vals_to_subs)<=0:
            # If the non-vero eigenavalues associated with the Jacobian matrix  
            # of non-DFE_1 are <= 0, the equil_pops is appended with non-DFE_1, 
            # as it is locally stable.
            equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
        
        # If all the values for non-DFE_2 are >= 0, evaluate eigenavalues 
        # associated with the Jacobian matrix of non-DFE_2.
        # Need to be able to subsititute in equilibria formula, as well as the param values:
        vals_to_subs = copy.deepcopy(param_values)
        vals_to_subs['S'] =  S_non_DFE_2
        vals_to_subs['I_W'] = I_W_non_DFE_2
        vals_to_subs['I_Z'] = I_Z_non_DFE_2
        if eigs_J_non_DFE_2[0].subs(vals_to_subs) <=0 and eigs_J_non_DFE_2[1].subs(vals_to_subs)<=0:
            # If the non-vero eigenavalues associated with the Jacobian matrix  
            # of non-DFE_2 are <= 0, the equil_pops is appended with non-DFE_2, 
            # as it is locally stable.
            equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    return(equil_pops)
    
    
    

def replace_inf_end_equil(param_values):
    '''
    Calculates the non-disease free equilibria for the replacement infection model ,
    as out lined in Spicknall et al 2013. Reuturning any non-DFE that are 
    that are biologically reasonable and locally stable.
    '''
    # Note if both strains have R0 values <1 there is no point in using this function.
    
    equil_pops = []
    # Setup Equilibrium populations as an empty list. If any of the endemic 
    # equilibria are found to be biologically feasible and stable they are appended to the list.
    
    #Basic Reproductive numbers from Spicknall 2013
    R0_W = beta_W/(gamma*(1-epsilon)+gamma_T*epsilon)
    R0_Z = beta_Z/gamma
    
    #Formula for the non-Disease Free Equilibria has been worked out in a jupyter notebook
    #and pasted here. NOTE with the first non-DFE the I_W population = N-S (I_Z=0),
    # and with the second non-DFE the I_Z population = N-S (I_W=0).
    # Therefore, only the formule for S are needed (see jupyter notebook).
    S_non_DFE_1 = N*(-epsilon*gamma + epsilon*gamma_T + gamma)/beta_W
    S_non_DFE_1 = S_non_DFE_1.subs(param_values)
    I_W_non_DFE_1 = param_values['N']-S_non_DFE_1
    I_Z_non_DFE_1 = 0
    S_non_DFE_2 = N*gamma/beta_Z
    S_non_DFE_2 = S_non_DFE_2.subs(param_values)
    I_W_non_DFE_2 = 0
    I_Z_non_DFE_2 = param_values['N']-S_non_DFE_2
    
    #The eigen values associated with jacobian matricies for the non-DFE have been 
    #worked out in a jupyter notebook and pasted here.
    eigs_J_non_DFE_1 = [-I_W*beta_W/N,
                        -(I_W*beta_W**2 - I_W*beta_W*beta_Z + N*beta_W*gamma + N*beta_Z*epsilon*gamma - N*beta_Z*epsilon*gamma_T - N*beta_Z*gamma)/(N*beta_W),
                        0]
    eigs_J_non_DFE_2 = [-I_Z*beta_Z/N,
                        (I_Z*beta_W*beta_Z - I_Z*beta_Z**2 + N*beta_W*gamma + N*beta_Z*epsilon*gamma - N*beta_Z*epsilon*gamma_T - N*beta_Z*gamma)/(N*beta_Z),
                        0]
    
    if math.isnan(S_non_DFE_1) or math.isnan(I_W_non_DFE_1) or math.isnan(I_Z_non_DFE_1) or math.isinf(S_non_DFE_1) or math.isinf(I_W_non_DFE_1) or math.isinf(I_Z_non_DFE_1):
        DFE_1_feasable = False
    else:
        DFE_1_feasable = True
    
    if DFE_1_feasable and round_sf(S_non_DFE_1,2) >= 0 and round_sf(I_W_non_DFE_1,2) >= 0 and round_sf(I_Z_non_DFE_1,2) >= 0:
        DFE_1_feasable = True
    else:
        DFE_1_feasable = False
    
    if math.isnan(S_non_DFE_2) or math.isnan(I_W_non_DFE_2) or math.isnan(I_Z_non_DFE_2) or math.isinf(S_non_DFE_2) or math.isinf(I_W_non_DFE_2) or math.isinf(I_Z_non_DFE_2):
        DFE_2_feasable = False
    else:
        DFE_2_feasable = True
    
    if DFE_2_feasable and round_sf(S_non_DFE_2,2) >= 0 and round_sf(I_W_non_DFE_2,2) >= 0 and round_sf(I_Z_non_DFE_2,2) >= 0:
        DFE_2_feasable = True
    else:
        DFE_2_feasable = False
    
    #Determine if R0_W is >=1, R0_Z<1 and the first DFE is biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)<1 and DFE_1_feasable:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_Z is >=1, R0_W<1 and the second DFE is biologically feasable.
    if R0_Z.subs(param_values)>=1 and R0_W.subs(param_values)<1 and DFE_2_feasable:   
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the first DFE is biologically feasable but the second is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable == False:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the second DFE is biologically feasable but the first is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable  == False and DFE_2_feasable:
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
        
    #Determine if R0_W is >=1, R0_Z>=1 and both DFEs are biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable:   
        # If all the values for non-DFE_1 are >= 0, evaluate eigenavalues 
        # associated with the Jacobian matrix of non-DFE_1.
        # Need to be able to subsititute in equilibria formula, as well as the param values:
        vals_to_subs = copy.deepcopy(param_values)
        vals_to_subs['S'] = S_non_DFE_1
        vals_to_subs['I_W'] = I_W_non_DFE_1
        vals_to_subs['I_Z'] = I_Z_non_DFE_1
        if eigs_J_non_DFE_1[0].subs(vals_to_subs) <=0 and eigs_J_non_DFE_1[1].subs(vals_to_subs)<=0:
            # If the non-vero eigenavalues associated with the Jacobian matrix  
            # of non-DFE_1 are <= 0, the equil_pops is appended with non-DFE_1, 
            # as it is locally stable.
            equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
        
        # If all the values for non-DFE_2 are >= 0, evaluate eigenavalues 
        # associated with the Jacobian matrix of non-DFE_2.
        # Need to be able to subsititute in equilibria formula, as well as the param values:
        vals_to_subs = copy.deepcopy(param_values)
        vals_to_subs['S'] =  S_non_DFE_2
        vals_to_subs['I_W'] = I_W_non_DFE_2
        vals_to_subs['I_Z'] = I_Z_non_DFE_2
        if eigs_J_non_DFE_2[0].subs(vals_to_subs) <=0 and eigs_J_non_DFE_2[1].subs(vals_to_subs)<=0:
            # If the non-vero eigenavalues associated with the Jacobian matrix  
            # of non-DFE_2 are <= 0, the equil_pops is appended with non-DFE_2, 
            # as it is locally stable.
            equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    return(equil_pops)
    
def possible_bi_directional_ee(param_values):
    '''
    Calculates the two possible endemic equilibria for the bi-directional conversion model,
    as out lined in Spicknall et al 2013.
    '''
    S_non_DFE_1 = N*(-beta_W*epsilon*phi + beta_W*gamma + beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*beta_Z)
    S_non_DFE_1 = S_non_DFE_1.subs(param_values)
    not_S_non_DFE_1 = param_values['N']-S_non_DFE_1
    I_W_non_DFE_1 = -I_Z*(beta_W*epsilon*phi - beta_W*gamma - beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*epsilon*rho)
    prop_I_W_non_DFE_1 = I_W_non_DFE_1.subs(param_values)/(1+I_W_non_DFE_1.subs(param_values))
    EE_1 = {'S': S_non_DFE_1.subs(param_values),
            'I_W': not_S_non_DFE_1*prop_I_W_non_DFE_1,
            'I_Z': not_S_non_DFE_1*(1-prop_I_W_non_DFE_1)}

    S_non_DFE_2 = N*(-beta_W*epsilon*phi + beta_W*gamma + beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*beta_Z)
    S_non_DFE_2 = S_non_DFE_2.subs(param_values)
    not_S_non_DFE_2 = param_values['N']-S_non_DFE_2
    I_W_non_DFE_2 = -I_Z*(beta_W*epsilon*phi - beta_W*gamma - beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*epsilon*rho)
    prop_I_W_non_DFE_2 = I_W_non_DFE_2.subs(param_values)/(1+I_W_non_DFE_2.subs(param_values))
    EE_2 = {'S': S_non_DFE_2.subs(param_values),
            'I_W': not_S_non_DFE_2*prop_I_W_non_DFE_2,
            'I_Z': not_S_non_DFE_2*(1-prop_I_W_non_DFE_2)}
    return {'Endemic Equilibria 1': EE_1, 'Endemic Equilibria 2': EE_2}

def bi_directional_end_equil(param_values):
    '''
    Calculates the non-disease free equilibria for the uni and bi-directional convertion model ,
    as out lined in Spicknall et al 2013. Returning any non-DFE that are 
    that are biologically reasonable and locally stable.
    '''
    # Note if both strains have R0 values <1 there is no point in using this function.
    
    equil_pops = []
    # Setup Equilibrium populations as an empty list. If any of the endemic 
    # equilibria are found to be biologically feasible they are appended to the list.
    
    #Basic Reproductive numbers from Spicknall 2013
    R0_W = 2*beta_W*beta_Z/(beta_W*gamma+beta_Z*gamma-beta_Z*epsilon*gamma+beta_Z*epsilon*gamma_T+beta_W*phi-beta_W*epsilon*phi+beta_Z*epsilon*rho-
                            ((beta_W*(gamma+phi-epsilon*phi)+beta_Z*(gamma-epsilon*gamma+epsilon*(gamma_T+rho)))**2+4*beta_W*beta_Z*((-1+epsilon)*gamma**2+
                                                                                                                                    (-1+epsilon)*epsilon*gamma_T*phi-gamma*(phi+phi*epsilon**2+epsilon*(gamma_T-2*phi+rho))))**0.5)
    R0_Z = 2*beta_W*beta_Z/(beta_W*gamma+beta_Z*gamma-beta_Z*epsilon*gamma+beta_Z*epsilon*gamma_T+beta_W*phi-beta_W*epsilon*phi+beta_Z*epsilon*rho+
                            ((beta_W*(gamma+phi-epsilon*phi)+beta_Z*(gamma-epsilon*gamma+epsilon*(gamma_T+rho)))**2+4*beta_W*beta_Z*((-1+epsilon)*gamma**2+
                                                                                                                                    (-1+epsilon)*epsilon*gamma_T*phi-gamma*(phi+phi*epsilon**2+epsilon*(gamma_T-2*phi+rho))))**0.5)
    
    #Formula for the non-Disease Free Equilibria has been worked out in a jupyter notebook
    #and pasted here.
    non_DFE_1 = {
            'S':N*(-beta_W*epsilon*phi + beta_W*gamma + beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*beta_Z)
            ,
            'I_W':-I_Z*(beta_W*epsilon*phi - beta_W*gamma - beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*epsilon*rho)
            ,
            'I_Z':I_Z
            }
    
    # Need to track proportion so a value for I_Z=1 needs to substituted along with params values.
    vals_to_subs = copy.deepcopy(param_values)
    vals_to_subs['I_Z'] = 1
    S_non_DFE_1 = non_DFE_1['S'].subs(param_values)
    not_S_non_DFE_1 = param_values['N']-S_non_DFE_1
    prop_I_W_non_DFE_1 = non_DFE_1['I_W'].subs(vals_to_subs)/(1+non_DFE_1['I_W'].subs(vals_to_subs))
    if math.isnan(prop_I_W_non_DFE_1):
        prop_I_W_non_DFE_1 =1
    I_W_non_DFE_1 = (not_S_non_DFE_1*prop_I_W_non_DFE_1)
    I_Z_non_DFE_1 = (not_S_non_DFE_1*(1-prop_I_W_non_DFE_1))
    
    non_DFE_2 = {
            'S':N*(-beta_W*epsilon*phi + beta_W*gamma + beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*beta_Z)
            ,
            'I_W':-I_Z*(beta_W*epsilon*phi - beta_W*gamma - beta_W*phi - beta_Z*epsilon*gamma + beta_Z*epsilon*gamma_T + beta_Z*epsilon*rho + beta_Z*gamma + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))/(2*beta_W*epsilon*rho)
            ,
            'I_Z':I_Z
            }
    
    S_non_DFE_2 = non_DFE_2['S'].subs(param_values)
    not_S_non_DFE_2 = param_values['N']-S_non_DFE_2
    prop_I_W_non_DFE_2 = non_DFE_2['I_W'].subs(vals_to_subs)/(1+non_DFE_2['I_W'].subs(vals_to_subs))
    if math.isnan(prop_I_W_non_DFE_2):
        prop_I_W_non_DFE_2 =1
    I_W_non_DFE_2 = (not_S_non_DFE_2*prop_I_W_non_DFE_2)
    I_Z_non_DFE_2 = (not_S_non_DFE_2*(1-prop_I_W_non_DFE_2))
    
        
    #The Eigen values associated with jacobian matricies for the non-DFE have been 
    #derived in a jupyter notebook and pasted here.
    eigs_J_non_DFE_1 = [I_Z*beta_W*phi/(4*N*rho) - I_Z*beta_W*gamma/(4*N*epsilon*rho) - I_Z*beta_W*phi/(4*N*epsilon*rho) - I_Z*beta_Z*gamma/(4*N*rho) + I_Z*beta_Z*gamma_T/(4*N*rho) - I_Z*beta_Z/(4*N) + I_Z*beta_Z*gamma/(4*N*epsilon*rho) - I_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*N*epsilon*rho) - beta_W*epsilon*phi/(4*beta_Z) + beta_W*gamma/(4*beta_Z) + beta_W*phi/(4*beta_Z) + epsilon*gamma/4 - epsilon*gamma_T/4 + epsilon*phi/4 - epsilon*rho/4 - gamma/2 - phi/4 - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_Z) - beta_Z*epsilon*gamma/(4*beta_W) + beta_Z*epsilon*gamma_T/(4*beta_W) + beta_Z*epsilon*rho/(4*beta_W) + beta_Z*gamma/(4*beta_W) - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_W) - sqrt(-8*I_Z*N*beta_W*beta_Z**2*epsilon*rho*(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi - beta_W*epsilon*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + 2*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_W*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_W*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2 + beta_Z*epsilon*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*epsilon*gamma_T*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)) + (-I_Z*beta_W**2*beta_Z*epsilon*phi + I_Z*beta_W**2*beta_Z*gamma + I_Z*beta_W**2*beta_Z*phi + I_Z*beta_W*beta_Z**2*epsilon*gamma - I_Z*beta_W*beta_Z**2*epsilon*gamma_T + I_Z*beta_W*beta_Z**2*epsilon*rho - I_Z*beta_W*beta_Z**2*gamma + I_Z*beta_W*beta_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + N*beta_W**2*epsilon**2*phi*rho - N*beta_W**2*epsilon*gamma*rho - N*beta_W**2*epsilon*phi*rho - N*beta_W*beta_Z*epsilon**2*gamma*rho + N*beta_W*beta_Z*epsilon**2*gamma_T*rho - N*beta_W*beta_Z*epsilon**2*phi*rho + N*beta_W*beta_Z*epsilon**2*rho**2 + 2*N*beta_W*beta_Z*epsilon*gamma*rho + N*beta_W*beta_Z*epsilon*phi*rho + N*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + N*beta_Z**2*epsilon**2*gamma*rho - N*beta_Z**2*epsilon**2*gamma_T*rho - N*beta_Z**2*epsilon**2*rho**2 - N*beta_Z**2*epsilon*gamma*rho + N*beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))**2)/(4*N*beta_W*beta_Z*epsilon*rho),
                        I_Z*beta_W*phi/(4*N*rho) - I_Z*beta_W*gamma/(4*N*epsilon*rho) - I_Z*beta_W*phi/(4*N*epsilon*rho) - I_Z*beta_Z*gamma/(4*N*rho) + I_Z*beta_Z*gamma_T/(4*N*rho) - I_Z*beta_Z/(4*N) + I_Z*beta_Z*gamma/(4*N*epsilon*rho) - I_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*N*epsilon*rho) - beta_W*epsilon*phi/(4*beta_Z) + beta_W*gamma/(4*beta_Z) + beta_W*phi/(4*beta_Z) + epsilon*gamma/4 - epsilon*gamma_T/4 + epsilon*phi/4 - epsilon*rho/4 - gamma/2 - phi/4 - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_Z) - beta_Z*epsilon*gamma/(4*beta_W) + beta_Z*epsilon*gamma_T/(4*beta_W) + beta_Z*epsilon*rho/(4*beta_W) + beta_Z*gamma/(4*beta_W) - sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_W) + sqrt(-8*I_Z*N*beta_W*beta_Z**2*epsilon*rho*(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi - beta_W*epsilon*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + 2*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_W*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_W*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2 + beta_Z*epsilon*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*epsilon*gamma_T*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_Z*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)) + (-I_Z*beta_W**2*beta_Z*epsilon*phi + I_Z*beta_W**2*beta_Z*gamma + I_Z*beta_W**2*beta_Z*phi + I_Z*beta_W*beta_Z**2*epsilon*gamma - I_Z*beta_W*beta_Z**2*epsilon*gamma_T + I_Z*beta_W*beta_Z**2*epsilon*rho - I_Z*beta_W*beta_Z**2*gamma + I_Z*beta_W*beta_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + N*beta_W**2*epsilon**2*phi*rho - N*beta_W**2*epsilon*gamma*rho - N*beta_W**2*epsilon*phi*rho - N*beta_W*beta_Z*epsilon**2*gamma*rho + N*beta_W*beta_Z*epsilon**2*gamma_T*rho - N*beta_W*beta_Z*epsilon**2*phi*rho + N*beta_W*beta_Z*epsilon**2*rho**2 + 2*N*beta_W*beta_Z*epsilon*gamma*rho + N*beta_W*beta_Z*epsilon*phi*rho + N*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + N*beta_Z**2*epsilon**2*gamma*rho - N*beta_Z**2*epsilon**2*gamma_T*rho - N*beta_Z**2*epsilon**2*rho**2 - N*beta_Z**2*epsilon*gamma*rho + N*beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))**2)/(4*N*beta_W*beta_Z*epsilon*rho),
                        0]
    eigs_J_non_DFE_2 = [I_Z*beta_W*phi/(4*N*rho) - I_Z*beta_W*gamma/(4*N*epsilon*rho) - I_Z*beta_W*phi/(4*N*epsilon*rho) - I_Z*beta_Z*gamma/(4*N*rho) + I_Z*beta_Z*gamma_T/(4*N*rho) - I_Z*beta_Z/(4*N) + I_Z*beta_Z*gamma/(4*N*epsilon*rho) + I_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*N*epsilon*rho) - beta_W*epsilon*phi/(4*beta_Z) + beta_W*gamma/(4*beta_Z) + beta_W*phi/(4*beta_Z) + epsilon*gamma/4 - epsilon*gamma_T/4 + epsilon*phi/4 - epsilon*rho/4 - gamma/2 - phi/4 + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_Z) - beta_Z*epsilon*gamma/(4*beta_W) + beta_Z*epsilon*gamma_T/(4*beta_W) + beta_Z*epsilon*rho/(4*beta_W) + beta_Z*gamma/(4*beta_W) + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_W) - sqrt(-8*I_Z*N*beta_W*beta_Z**2*epsilon*rho*(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_W*epsilon*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - 2*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_W*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_W*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2 - beta_Z*epsilon*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*epsilon*gamma_T*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)) + (I_Z*beta_W**2*beta_Z*epsilon*phi - I_Z*beta_W**2*beta_Z*gamma - I_Z*beta_W**2*beta_Z*phi - I_Z*beta_W*beta_Z**2*epsilon*gamma + I_Z*beta_W*beta_Z**2*epsilon*gamma_T - I_Z*beta_W*beta_Z**2*epsilon*rho + I_Z*beta_W*beta_Z**2*gamma + I_Z*beta_W*beta_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - N*beta_W**2*epsilon**2*phi*rho + N*beta_W**2*epsilon*gamma*rho + N*beta_W**2*epsilon*phi*rho + N*beta_W*beta_Z*epsilon**2*gamma*rho - N*beta_W*beta_Z*epsilon**2*gamma_T*rho + N*beta_W*beta_Z*epsilon**2*phi*rho - N*beta_W*beta_Z*epsilon**2*rho**2 - 2*N*beta_W*beta_Z*epsilon*gamma*rho - N*beta_W*beta_Z*epsilon*phi*rho + N*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - N*beta_Z**2*epsilon**2*gamma*rho + N*beta_Z**2*epsilon**2*gamma_T*rho + N*beta_Z**2*epsilon**2*rho**2 + N*beta_Z**2*epsilon*gamma*rho + N*beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))**2)/(4*N*beta_W*beta_Z*epsilon*rho),
                        I_Z*beta_W*phi/(4*N*rho) - I_Z*beta_W*gamma/(4*N*epsilon*rho) - I_Z*beta_W*phi/(4*N*epsilon*rho) - I_Z*beta_Z*gamma/(4*N*rho) + I_Z*beta_Z*gamma_T/(4*N*rho) - I_Z*beta_Z/(4*N) + I_Z*beta_Z*gamma/(4*N*epsilon*rho) + I_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*N*epsilon*rho) - beta_W*epsilon*phi/(4*beta_Z) + beta_W*gamma/(4*beta_Z) + beta_W*phi/(4*beta_Z) + epsilon*gamma/4 - epsilon*gamma_T/4 + epsilon*phi/4 - epsilon*rho/4 - gamma/2 - phi/4 + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_Z) - beta_Z*epsilon*gamma/(4*beta_W) + beta_Z*epsilon*gamma_T/(4*beta_W) + beta_Z*epsilon*rho/(4*beta_W) + beta_Z*gamma/(4*beta_W) + sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)/(4*beta_W) + sqrt(-8*I_Z*N*beta_W*beta_Z**2*epsilon*rho*(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_W*epsilon*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - 2*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_W*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - beta_W*phi*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2 - beta_Z*epsilon*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*epsilon*gamma_T*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) + beta_Z*gamma*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2)) + (I_Z*beta_W**2*beta_Z*epsilon*phi - I_Z*beta_W**2*beta_Z*gamma - I_Z*beta_W**2*beta_Z*phi - I_Z*beta_W*beta_Z**2*epsilon*gamma + I_Z*beta_W*beta_Z**2*epsilon*gamma_T - I_Z*beta_W*beta_Z**2*epsilon*rho + I_Z*beta_W*beta_Z**2*gamma + I_Z*beta_W*beta_Z*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - N*beta_W**2*epsilon**2*phi*rho + N*beta_W**2*epsilon*gamma*rho + N*beta_W**2*epsilon*phi*rho + N*beta_W*beta_Z*epsilon**2*gamma*rho - N*beta_W*beta_Z*epsilon**2*gamma_T*rho + N*beta_W*beta_Z*epsilon**2*phi*rho - N*beta_W*beta_Z*epsilon**2*rho**2 - 2*N*beta_W*beta_Z*epsilon*gamma*rho - N*beta_W*beta_Z*epsilon*phi*rho + N*beta_W*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2) - N*beta_Z**2*epsilon**2*gamma*rho + N*beta_Z**2*epsilon**2*gamma_T*rho + N*beta_Z**2*epsilon**2*rho**2 + N*beta_Z**2*epsilon*gamma*rho + N*beta_Z*epsilon*rho*sqrt(beta_W**2*epsilon**2*phi**2 - 2*beta_W**2*epsilon*gamma*phi - 2*beta_W**2*epsilon*phi**2 + beta_W**2*gamma**2 + 2*beta_W**2*gamma*phi + beta_W**2*phi**2 - 2*beta_W*beta_Z*epsilon**2*gamma*phi + 2*beta_W*beta_Z*epsilon**2*gamma_T*phi - 2*beta_W*beta_Z*epsilon**2*phi*rho + 2*beta_W*beta_Z*epsilon*gamma**2 - 2*beta_W*beta_Z*epsilon*gamma*gamma_T + 4*beta_W*beta_Z*epsilon*gamma*phi - 2*beta_W*beta_Z*epsilon*gamma*rho - 2*beta_W*beta_Z*epsilon*gamma_T*phi + 2*beta_W*beta_Z*epsilon*phi*rho - 2*beta_W*beta_Z*gamma**2 - 2*beta_W*beta_Z*gamma*phi + beta_Z**2*epsilon**2*gamma**2 - 2*beta_Z**2*epsilon**2*gamma*gamma_T - 2*beta_Z**2*epsilon**2*gamma*rho + beta_Z**2*epsilon**2*gamma_T**2 + 2*beta_Z**2*epsilon**2*gamma_T*rho + beta_Z**2*epsilon**2*rho**2 - 2*beta_Z**2*epsilon*gamma**2 + 2*beta_Z**2*epsilon*gamma*gamma_T + 2*beta_Z**2*epsilon*gamma*rho + beta_Z**2*gamma**2))**2)/(4*N*beta_W*beta_Z*epsilon*rho),
                        0]
    
    
    if math.isnan(S_non_DFE_1) or math.isnan(I_W_non_DFE_1) or math.isnan(I_Z_non_DFE_1) or math.isinf(S_non_DFE_1) or math.isinf(I_W_non_DFE_1) or math.isinf(I_Z_non_DFE_1):
        DFE_1_feasable = False
    else:
        DFE_1_feasable = True
    
    if DFE_1_feasable and round_sf(S_non_DFE_1,2) >= 0 and round_sf(I_W_non_DFE_1,2) >= 0 and round_sf(I_Z_non_DFE_1,2) >= 0:
        DFE_1_feasable = True
    else:
        DFE_1_feasable = False
    
    if math.isnan(S_non_DFE_2) or math.isnan(I_W_non_DFE_2) or math.isnan(I_Z_non_DFE_2) or math.isinf(S_non_DFE_2) or math.isinf(I_W_non_DFE_2) or math.isinf(I_Z_non_DFE_2):
        DFE_2_feasable = False
    else:
        DFE_2_feasable = True
    
    if DFE_2_feasable and round_sf(S_non_DFE_2,2) >= 0 and round_sf(I_W_non_DFE_2,2) >= 0 and round_sf(I_Z_non_DFE_2,2) >= 0:
        DFE_2_feasable = True
    else:
        DFE_2_feasable = False
    
    #Determine if R0_W is >=1, R0_Z<1 and the first DFE is biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)<1 and DFE_1_feasable:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_Z is >=1, R0_W<1 and the second DFE is biologically feasable.
    if R0_Z.subs(param_values)>=1 and R0_W.subs(param_values)<1 and DFE_2_feasable:   
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the first DFE is biologically feasable but the second is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable == False:
        equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
    
    #Determine if R0_W is >=1, R0_Z>=1 and the second DFE is biologically feasable but the first is not.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable  == False and DFE_2_feasable:
        equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
        
    #Determine if R0_W is >=1, R0_Z>=1 and both DFEs are biologically feasable.
    if R0_W.subs(param_values)>=1 and R0_Z.subs(param_values)>=1 and DFE_1_feasable and DFE_2_feasable:
        # There are several denominators in eigen values associated with both DFEs (see jupyter notebook) that can be 0 if rho or epsilon = 0.
        # This means sympy's substitution function can produce undefined values (NaN values).
        # Therefore as Spicknall et al (2013) found that the largest of the R0 values was the R0 that described the system.
        # We can choose between the DFE's by using R0s.
        if param_values['epsilon']==0 or param_values['rho'] == 0:
            if R0_W.subs(param_values)> R0_Z.subs(param_values):
                equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
            else:
                equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
        else:
            # Need to be able to subsititute in equilibria formula, as well as the param values:
            vals_to_subs = copy.deepcopy(param_values)
            vals_to_subs['S'] = S_non_DFE_1
            vals_to_subs['I_W'] = I_W_non_DFE_1
            vals_to_subs['I_Z'] = I_Z_non_DFE_1
            # Both DFE associated eigen values can be complex and we are only 
            # interested in the real part so:
            eig_1_J_non_DFE_1 = sympy.re(eigs_J_non_DFE_1[0].subs(vals_to_subs))
            eig_2_J_non_DFE_1 = sympy.re(eigs_J_non_DFE_1[1].subs(vals_to_subs))
            if eig_1_J_non_DFE_1 <=0 and eig_2_J_non_DFE_1 <=0:
                # If the non-vero eigenavalues associated with the Jacobian matrix  
                # of non-DFE_1 are <= 0, the equil_pops is appended with non-DFE_1, 
                # as it is locally stable.
                equil_pops.append({'S':S_non_DFE_1,'I_W':I_W_non_DFE_1,'I_Z':I_Z_non_DFE_1})
            
            # If all the values for non-DFE_2 are >= 0, evaluate eigenavalues 
            # associated with the Jacobian matrix of non-DFE_2.
            # Need to be able to subsititute in equilibria formula, as well as the param values:
            vals_to_subs = copy.deepcopy(param_values)
            vals_to_subs['S'] =  S_non_DFE_2
            vals_to_subs['I_W'] = I_W_non_DFE_2
            vals_to_subs['I_Z'] = I_Z_non_DFE_2
            # Both DFE associated eigen values can be complex and we are only 
            # interested in the real part so:
            eig_1_J_non_DFE_2 = sympy.re(eigs_J_non_DFE_2[0].subs(vals_to_subs))
            eig_2_J_non_DFE_2 = sympy.re(eigs_J_non_DFE_2[1].subs(vals_to_subs))
            if eig_1_J_non_DFE_2 <=0 and eig_2_J_non_DFE_2 <=0:
                # If the non-vero eigenavalues associated with the Jacobian matrix  
                # of non-DFE_2 are <= 0, the equil_pops is appended with non-DFE_2, 
                # as it is locally stable.
                equil_pops.append({'S':S_non_DFE_2,'I_W':I_W_non_DFE_2,'I_Z':I_Z_non_DFE_2})
    
    return(equil_pops)