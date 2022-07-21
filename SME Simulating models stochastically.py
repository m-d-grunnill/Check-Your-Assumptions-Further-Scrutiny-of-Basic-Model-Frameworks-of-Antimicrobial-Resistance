"""
Creation:
    Author: Martin Grunnill
    Date: 21/07/2022
Description: Code for running models stochastically and producing figures comparing the outcomes from the different
             scenarios and deterministc simulation.
             Note you need to run SMC Simulating models deterministically.ipynb in order to have the deterministic
             data for producing comparison figures.
"""
import os
import psutil
import datetime
import numpy as np
from timeit import default_timer as timer
from dask.distributed import Client
from tqdm import tqdm, trange
import SMB_Code_models_equilibria as models
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#%% Functions used in seting up and running stochastic simulation scenarios.

def exclus_n_replace_inf_I_W_equilibrium(param_values):
    '''
    Calculates the equilibia for where the sensitive strain  exists for the
    exclusive and replacement infection models,
    as out lined in Spicknall et al 2013.     '''

    # Formula for the non-Disease Free Equilibria has been worked out in a jupyter notebook
    # and pasted here. NOTE for this model only the susceptible variables equations is needed.
    S_non_DFE_1 = (-param_values['epsilon'] * param_values['gamma'] + param_values['epsilon'] * param_values['gamma_T']
                   + param_values['gamma']) / (param_values['beta_W'] / param_values['N'])
    S_non_DFE_1 = int(round(S_non_DFE_1))
    I_W_non_DFE_1 = int(round(param_values['N'] - S_non_DFE_1))
    I_Z_non_DFE_1 = 0
    equil_props = {'S': S_non_DFE_1, 'I_W': I_W_non_DFE_1, 'I_Z': I_Z_non_DFE_1}
    return equil_props

def full_co_inf_W_equilibrium(param_values):
    '''
    Calculates the endemic for the existance sensitive strain  only in
    full coinfection model.
    '''
    EE_S = (param_values['N'] *
            (-param_values['epsilon'] * param_values['gamma'] + param_values['epsilon'] * param_values['gamma_T'] + param_values['gamma']) /
            param_values['beta_W'])
    EE_I_W = (param_values['N'] *
              (-param_values['epsilon'] * param_values['gamma'] + param_values['epsilon'] * param_values['gamma_T'] + param_values['gamma']) *
              (param_values['beta_W'] + param_values['epsilon'] * param_values['gamma'] - param_values['epsilon'] * param_values['gamma_T'] - param_values['gamma']) / param_values['beta_W'] ** 2)
    EE_I_WW = (param_values['N'] *
               (param_values['beta_W'] ** 2 - 2 * param_values['beta_W'] *
                (-param_values['epsilon'] * param_values['gamma'] + param_values['epsilon'] * param_values['gamma_T'] + param_values['gamma']) +
                param_values['epsilon'] ** 2 * param_values['gamma'] ** 2 - 2 * param_values['epsilon'] ** 2 *
                param_values['gamma'] * param_values['gamma_T'] + param_values['epsilon'] ** 2 *
                param_values['gamma_T'] ** 2 - 2 * param_values['epsilon'] * param_values['gamma'] ** 2 +
                2 * param_values['epsilon'] * param_values['gamma'] * param_values['gamma_T'] + param_values['gamma'] ** 2)
               / param_values['beta_W'] ** 2)
    EE_S = int(round(EE_S))
    EE_I_W = int(round(EE_I_W))
    EE_I_WW = int(round(EE_I_WW))
    return {'S': EE_S, 'I_W': EE_I_W, 'I_WW': EE_I_WW, 'I_Z': 0, 'I_ZZ': 0, 'I_WZ': 0}


def round_sf(number, significant):
    return round(number, significant - len(str(number)))

def run_model_stoch_and_det(model, model_name, save_dir, time_range, print_msg, iterations=100, stoch_chunks=10):
    stoch_data_name = save_dir + '/' + 'Stochastic data'
    iterations_per_chunk = iterations / stoch_chunks
    if not iterations_per_chunk.is_integer():
        raise ValueError('iterations divided by stoch_chunks must produce no remainder.')
    iterations_per_chunk = int(iterations_per_chunk)
    if os.path.exists(save_dir + '/' + 'Coinfection stochastic data.npy'):
        os.rename(save_dir + '/' + 'Coinfection stochastic data.npy', stoch_data_name + '.npy')

    if not os.path.exists(stoch_data_name + '.npy'):
        now = datetime.datetime.now()
        print(model_name + ' ' + print_msg + ' stoch sims started at ' + now.strftime(
            "%Y-%m-%d %H:%M:%S") + '.')
        start_time = timer()
        results = []
        for i in trange(stoch_chunks, desc='Chunk of stochastic simulations'):
            simX, simT = model.simulate_jump(time_range, iteration=iterations_per_chunk, full_output=True, parallel=True)
            results += simX
        results = np.dstack(results)
        np.save(stoch_data_name, results)
        end_time = timer()
        time_of_run = end_time - start_time
        now = datetime.datetime.now()
        processing_time = 'Time of run =' + str(round(time_of_run, 2) / 60) + 'minutes.'
        with open(save_dir + '/run time.txt', 'w') as f:
            f.write(processing_time)
        print(model_name + ' ' + print_msg + ' stoch sims ended at ' + now.strftime("%Y-%m-%d %H:%M:%S") +
              '.\n' + processing_time)
    det_data_name = save_dir + '/' + 'Detereministic data'
    if not os.path.exists(det_data_name + '.npy'):
        det_run = model.integrate(time_range[1:])
        np.save(det_data_name, det_run)

def run_model_with_scenario_1(model, param_values, epsilon_range, model_name, directory, time_range):
    mod_directory = directory+'/'+model_name
    if not os.path.exists(mod_directory):
        os.mkdir(mod_directory)
    state_index_dict = {state.ID: index for index, state in enumerate(model.state_list)}
    for epsilon_value in tqdm(epsilon_range, desc=model_name + ' ith in epsilon_range'):
        epsilon_value_str = str(epsilon_value)
        param_values['epsilon'] = epsilon_value
        init_state = np.zeros(model.num_state)
        if model_name == 'Full coinfection model':
            AMS_equil = full_co_inf_W_equilibrium(param_values)
            init_state[state_index_dict['I_WW']] = AMS_equil['I_WW']
        else:
            AMS_equil = exclus_n_replace_inf_I_W_equilibrium(param_values)

        S_initial = AMS_equil['S']-int(1)
        if model_name == 'Single strain model':
            init_state[state_index_dict['S_W']] = S_initial
            init_state[state_index_dict['S_Z']] = S_initial
        else:
            init_state[state_index_dict['S']] = S_initial

        init_state[state_index_dict['I_W']] = AMS_equil['I_W']
        init_state[state_index_dict['I_Z']] = int(1)
        model.initial_values = (init_state, time_range[0])
        model.parameters = param_values
        save_dir = mod_directory + '/epsilon value = ' + epsilon_value_str
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        print_msg = 'Epsilon value ' + epsilon_value_str
        run_model_stoch_and_det(model, model_name, save_dir, time_range, print_msg)

def run_model_with_scenario_2(model, param_values, epsilon_range, model_name, directory, time_range):
    mod_directory = directory+'/'+model_name
    if not os.path.exists(mod_directory):
        os.mkdir(mod_directory)
    state_index_dict = {state.ID: index for index, state in enumerate(model.state_list)}
    init_state = np.zeros(model.num_state)
    if model_name == 'Single strain model':
        init_state[state_index_dict['S_W']] = int(param_values['N'] - 2)
        init_state[state_index_dict['S_Z']] = int(param_values['N'] - 2)
    else:
        init_state[state_index_dict['S']] = int(param_values['N']-2)
    init_state[state_index_dict['I_W']] = int(1)
    init_state[state_index_dict['I_Z']] = int(1)
    model.initial_values = (init_state ,  time_range[0])
    for epsilon_value in tqdm(epsilon_range, desc=model_name + ' ith in epsilon_range'):
        epsilon_value_str = str(epsilon_value)
        param_values['epsilon'] = epsilon_value
        model.parameters = param_values
        save_dir = mod_directory + '/epsilon value = ' + epsilon_value_str
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        print_msg = 'Epsilon value ' + epsilon_value_str
        run_model_stoch_and_det(model, model_name, save_dir, time_range, print_msg)

#%% Functions used to produce extra lines on graphs

# Formula for the non-Disease Free Equilibria has been worked out in a jupyter notebook
# and pasted here. NOTE with the first non-DFE the I_W population = N-S (I_Z=0),
# and with the second non-DFE the I_Z population = N-S (I_W=0).
# Therefore, only the formule for S are needed (see jupyter notebook).
def possible_exclusive_replacement_ee(N, beta_W, beta_Z, gamma, gamma_T, epsilon_values_as_percent):
    '''
    Calculates the two two possible endemic equilibria for the Exclusive and Replacement infection models,
    as out lined in Spicknall et al 2013. Reuturning any non-DFE that are
    that are biologically reasonable and locally stable.
    '''
    epsilon_values = np.divide(epsilon_values_as_percent,100)
    def S_EE_1(N, beta_W, beta_Z, gamma, gamma_T, epsilon):
        return N*(-epsilon*gamma + epsilon * gamma_T + gamma)/beta_W

    S_EE_1 = np.vectorize(S_EE_1)
    S_non_DFE_1 = S_EE_1(N, beta_W, beta_Z, gamma, gamma_T, epsilon_values)
    EE_1 = {'% Receiving Treatment': epsilon_values_as_percent,
            'Susceptible': S_non_DFE_1,
            'All Sensitive Infections': np.subtract(N,S_non_DFE_1),
            'All Resistant Infections': 0}
    S_non_DFE_2 = np.repeat(N*gamma/beta_Z,len(epsilon_values))
    EE_2 = {'% Receiving Treatment': epsilon_values_as_percent,
            'Susceptible': S_non_DFE_2,
            'All Sensitive Infections': 0,
            'All Resistant Infections': np.subtract(N,S_non_DFE_2)}
    return pd.DataFrame(EE_1), pd.DataFrame(EE_2)


def possible_bi_directional_ee(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon_values_as_percent):
    '''
    Calculates the two possible endemic equilibria for the bi-directional conversion model,
    as out lined in Spicknall et al 2013.
    '''
    epsilon_values = np.divide(epsilon_values_as_percent, 100)

    def get_S_ee1(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon):
        return N * (
                    -beta_W * epsilon * phi + beta_W * gamma + beta_W * phi - beta_Z * epsilon * gamma + beta_Z * epsilon * gamma_T + beta_Z * epsilon * rho + beta_Z * gamma - np.sqrt(
                beta_W ** 2 * epsilon ** 2 * phi ** 2 - 2 * beta_W ** 2 * epsilon * gamma * phi - 2 * beta_W ** 2 * epsilon * phi ** 2 + beta_W ** 2 * gamma ** 2 + 2 * beta_W ** 2 * gamma * phi + beta_W ** 2 * phi ** 2 - 2 * beta_W * beta_Z * epsilon ** 2 * gamma * phi + 2 * beta_W * beta_Z * epsilon ** 2 * gamma_T * phi - 2 * beta_W * beta_Z * epsilon ** 2 * phi * rho + 2 * beta_W * beta_Z * epsilon * gamma ** 2 - 2 * beta_W * beta_Z * epsilon * gamma * gamma_T + 4 * beta_W * beta_Z * epsilon * gamma * phi - 2 * beta_W * beta_Z * epsilon * gamma * rho - 2 * beta_W * beta_Z * epsilon * gamma_T * phi + 2 * beta_W * beta_Z * epsilon * phi * rho - 2 * beta_W * beta_Z * gamma ** 2 - 2 * beta_W * beta_Z * gamma * phi + beta_Z ** 2 * epsilon ** 2 * gamma ** 2 - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * gamma_T - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * rho + beta_Z ** 2 * epsilon ** 2 * gamma_T ** 2 + 2 * beta_Z ** 2 * epsilon ** 2 * gamma_T * rho + beta_Z ** 2 * epsilon ** 2 * rho ** 2 - 2 * beta_Z ** 2 * epsilon * gamma ** 2 + 2 * beta_Z ** 2 * epsilon * gamma * gamma_T + 2 * beta_Z ** 2 * epsilon * gamma * rho + beta_Z ** 2 * gamma ** 2)) / (
                           2 * beta_W * beta_Z)

    get_S_ee1 = np.vectorize(get_S_ee1)

    def get_prop_I_W_ee1(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon):
        I_W_ee1_sym = -(
                    beta_W * epsilon * phi - beta_W * gamma - beta_W * phi - beta_Z * epsilon * gamma + beta_Z * epsilon * gamma_T + beta_Z * epsilon * rho + beta_Z * gamma - np.sqrt(
                beta_W ** 2 * epsilon ** 2 * phi ** 2 - 2 * beta_W ** 2 * epsilon * gamma * phi - 2 * beta_W ** 2 * epsilon * phi ** 2 + beta_W ** 2 * gamma ** 2 + 2 * beta_W ** 2 * gamma * phi + beta_W ** 2 * phi ** 2 - 2 * beta_W * beta_Z * epsilon ** 2 * gamma * phi + 2 * beta_W * beta_Z * epsilon ** 2 * gamma_T * phi - 2 * beta_W * beta_Z * epsilon ** 2 * phi * rho + 2 * beta_W * beta_Z * epsilon * gamma ** 2 - 2 * beta_W * beta_Z * epsilon * gamma * gamma_T + 4 * beta_W * beta_Z * epsilon * gamma * phi - 2 * beta_W * beta_Z * epsilon * gamma * rho - 2 * beta_W * beta_Z * epsilon * gamma_T * phi + 2 * beta_W * beta_Z * epsilon * phi * rho - 2 * beta_W * beta_Z * gamma ** 2 - 2 * beta_W * beta_Z * gamma * phi + beta_Z ** 2 * epsilon ** 2 * gamma ** 2 - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * gamma_T - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * rho + beta_Z ** 2 * epsilon ** 2 * gamma_T ** 2 + 2 * beta_Z ** 2 * epsilon ** 2 * gamma_T * rho + beta_Z ** 2 * epsilon ** 2 * rho ** 2 - 2 * beta_Z ** 2 * epsilon * gamma ** 2 + 2 * beta_Z ** 2 * epsilon * gamma * gamma_T + 2 * beta_Z ** 2 * epsilon * gamma * rho + beta_Z ** 2 * gamma ** 2)) / (
                                  2 * beta_W * epsilon * rho)
        return I_W_ee1_sym / (1 + I_W_ee1_sym)

    get_prop_I_W_ee1 = np.vectorize(get_prop_I_W_ee1)
    epsilon_values = np.divide(epsilon_values_as_percent, 100)

    S_ee1 = get_S_ee1(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon_values)
    prop_I_W_ee1 = get_prop_I_W_ee1(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon_values)
    not_S_ee1_sym = np.subtract(N, S_ee1)
    EE1 = {'% Receiving Treatment': epsilon_values_as_percent,
           'Susceptible': S_ee1,
           'All Sensitive Infections': np.multiply(not_S_ee1_sym, prop_I_W_ee1),
           'All Resistant Infections': np.multiply(not_S_ee1_sym, np.subtract(1, prop_I_W_ee1))}

    def get_S_ee2(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon):
        return N * (
                    -beta_W * epsilon * phi + beta_W * gamma + beta_W * phi - beta_Z * epsilon * gamma + beta_Z * epsilon * gamma_T + beta_Z * epsilon * rho + beta_Z * gamma + np.sqrt(
                beta_W ** 2 * epsilon ** 2 * phi ** 2 - 2 * beta_W ** 2 * epsilon * gamma * phi - 2 * beta_W ** 2 * epsilon * phi ** 2 + beta_W ** 2 * gamma ** 2 + 2 * beta_W ** 2 * gamma * phi + beta_W ** 2 * phi ** 2 - 2 * beta_W * beta_Z * epsilon ** 2 * gamma * phi + 2 * beta_W * beta_Z * epsilon ** 2 * gamma_T * phi - 2 * beta_W * beta_Z * epsilon ** 2 * phi * rho + 2 * beta_W * beta_Z * epsilon * gamma ** 2 - 2 * beta_W * beta_Z * epsilon * gamma * gamma_T + 4 * beta_W * beta_Z * epsilon * gamma * phi - 2 * beta_W * beta_Z * epsilon * gamma * rho - 2 * beta_W * beta_Z * epsilon * gamma_T * phi + 2 * beta_W * beta_Z * epsilon * phi * rho - 2 * beta_W * beta_Z * gamma ** 2 - 2 * beta_W * beta_Z * gamma * phi + beta_Z ** 2 * epsilon ** 2 * gamma ** 2 - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * gamma_T - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * rho + beta_Z ** 2 * epsilon ** 2 * gamma_T ** 2 + 2 * beta_Z ** 2 * epsilon ** 2 * gamma_T * rho + beta_Z ** 2 * epsilon ** 2 * rho ** 2 - 2 * beta_Z ** 2 * epsilon * gamma ** 2 + 2 * beta_Z ** 2 * epsilon * gamma * gamma_T + 2 * beta_Z ** 2 * epsilon * gamma * rho + beta_Z ** 2 * gamma ** 2)) / (
                           2 * beta_W * beta_Z)

    get_S_ee2 = np.vectorize(get_S_ee2)

    def get_prop_I_W_ee2(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon):
        I_W_ee2_sym = -(
                    beta_W * epsilon * phi - beta_W * gamma - beta_W * phi - beta_Z * epsilon * gamma + beta_Z * epsilon * gamma_T + beta_Z * epsilon * rho + beta_Z * gamma + np.sqrt(
                beta_W ** 2 * epsilon ** 2 * phi ** 2 - 2 * beta_W ** 2 * epsilon * gamma * phi - 2 * beta_W ** 2 * epsilon * phi ** 2 + beta_W ** 2 * gamma ** 2 + 2 * beta_W ** 2 * gamma * phi + beta_W ** 2 * phi ** 2 - 2 * beta_W * beta_Z * epsilon ** 2 * gamma * phi + 2 * beta_W * beta_Z * epsilon ** 2 * gamma_T * phi - 2 * beta_W * beta_Z * epsilon ** 2 * phi * rho + 2 * beta_W * beta_Z * epsilon * gamma ** 2 - 2 * beta_W * beta_Z * epsilon * gamma * gamma_T + 4 * beta_W * beta_Z * epsilon * gamma * phi - 2 * beta_W * beta_Z * epsilon * gamma * rho - 2 * beta_W * beta_Z * epsilon * gamma_T * phi + 2 * beta_W * beta_Z * epsilon * phi * rho - 2 * beta_W * beta_Z * gamma ** 2 - 2 * beta_W * beta_Z * gamma * phi + beta_Z ** 2 * epsilon ** 2 * gamma ** 2 - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * gamma_T - 2 * beta_Z ** 2 * epsilon ** 2 * gamma * rho + beta_Z ** 2 * epsilon ** 2 * gamma_T ** 2 + 2 * beta_Z ** 2 * epsilon ** 2 * gamma_T * rho + beta_Z ** 2 * epsilon ** 2 * rho ** 2 - 2 * beta_Z ** 2 * epsilon * gamma ** 2 + 2 * beta_Z ** 2 * epsilon * gamma * gamma_T + 2 * beta_Z ** 2 * epsilon * gamma * rho + beta_Z ** 2 * gamma ** 2)) / (
                                  2 * beta_W * epsilon * rho)
        return I_W_ee2_sym / (1 + I_W_ee2_sym)

    get_prop_I_W_ee2 = np.vectorize(get_prop_I_W_ee2)

    S_ee2 = get_S_ee2(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon_values)
    prop_I_W_ee2 = get_prop_I_W_ee2(N, beta_W, beta_Z, gamma, gamma_T, rho, phi, epsilon_values)
    not_S_ee2_sym = np.subtract(N, S_ee2)
    EE2 = {'% Receiving Treatment': epsilon_values_as_percent,
           'Susceptible': S_ee2,
           'All Sensitive Infections': np.multiply(not_S_ee2_sym, prop_I_W_ee2),
           'All Resistant Infections': np.multiply(not_S_ee2_sym, np.subtract(1, prop_I_W_ee2))}
    return pd.DataFrame(EE1).dropna(), pd.DataFrame(EE2).dropna()

if __name__ == '__main__':
    #%% Setting parrallel simulation and save directory
    # save directory
    data_directory = 'D:/AMR models'
    if not os.path.exists(data_directory):
        os.mkdir(data_directory)

    # setup the dask cluster
    number_of_workers = os.cpu_count() # alter this if you do not want all your computing resources used.
    memory = psutil.virtual_memory() # alter this if you do not want all your computing resources used.
    memory_to_use = memory.total
    client = Client(n_workers=number_of_workers, threads_per_worker=1, memory_limit=memory_to_use)


    #%% Running stochastic simulations
    ## Setup population, parameters , iterrations, runtime and epsilon values
    years=100
    N=int(1e3)
    days = 365*years
    time_range = np.arange(0, days+1, 1)


    #%%
    # Scenario 2
    scn_2_directory = data_directory + '/Scenario 2'
    if not os.path.exists(scn_2_directory):
        os.mkdir(scn_2_directory)
    epsilon_range = [i / 100 for i in range(0, 51)]

    # %%
    param_values = {'N': int(N), 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1, 'q': 0.5}
    model_name = 'Full coinfection model'
    model = models.AMR_full_co_inf_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)

    # %%
    model_name = 'Superinfection model'
    model = models.AMR_superinf_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)

    #%%
    param_values = {'N': int(N), 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1}
    model_name = 'Exclusive infection model'
    model = models.AMR_exclus_inf_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)
    
    #%%
    model_name = 'Replacement infection model'
    model = models.AMR_replace_inf_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)

    #%%
    model_name = 'Unidirectional conversion model'
    param_values.update({'rho': 0.5, 'phi': 0.0})
    model = models.AMR_bi_and_uni_conversion_inf_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)

    #%%
    model_name = 'Bidirectional conversion model'
    param_values.update({'phi':0.05})
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)

    #%%
    model_name = 'Single strain model'
    param_values = {'N': int(N), 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1}
    model = models.AMR_single_strain_intia()
    run_model_with_scenario_2(model, param_values, epsilon_range, model_name, scn_2_directory, time_range)



    #%%
    # Scenario 1
    scn_1_directory = data_directory + '/Scenario 1'
    if not os.path.exists(scn_1_directory):
        os.mkdir(scn_1_directory)

    epsilon_range = [i / 100 for i in range(33, -1, -1)]
    param_values = {'N': int(N), 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1}
    #%%
    model_name = 'Exclusive infection model'
    model = models.AMR_exclus_inf_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)
    #%%
    model_name = 'Replacement infection model'
    model = models.AMR_replace_inf_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)

    #%%
    model_name = 'Unidirectional conversion model'
    param_values.update({'rho': 0.5, 'phi': 0.0})
    model = models.AMR_bi_and_uni_conversion_inf_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)

    #%%
    model_name = 'Bidirectional conversion model'
    param_values.update({'phi':0.05})
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)

    #%%
    # single strain model
    model_name = 'Single strain model'
    param_values = {'N': int(N), 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1}
    model = models.AMR_single_strain_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)

    #%%
    model_name = 'Superinfection model'
    param_values= {'N':int(N),'beta_W': 0.04 , 'beta_Z': 0.015,'gamma': 0.01,'gamma_T': 0.1,'q':0.5}
    model = models.AMR_superinf_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)

    #%%
    model_name = 'Full coinfection model'
    model = models.AMR_full_co_inf_intia()
    run_model_with_scenario_1(model, param_values, epsilon_range, model_name, scn_1_directory, time_range)
    
    #%% Merging npy files.
    scenarios_epsilon_values = {'Scenario 1': [i / 100 for i in range(0, 34)],
                                'Scenario 2': [i / 100 for i in range(0, 51)]}

    models = ['Exclusive infection model', 'Replacement infection model',
              'Superinfection model', 'Full coinfection model',
              'Bidirectional conversion model', 'Unidirectional conversion model',
              'Single strain model']

    for scenario, epsilon_values in tqdm.tqdm(scenarios_epsilon_values.items(), desc='Scenario'):
        scenario_path = data_directory + '/' + scenario
        for model in tqdm.tqdm(models, desc='model'):
            deteministic_data = []
            stochastic_data = []
            model_path = scenario_path + '/' + model
            for epsilon_value in tqdm.tqdm(epsilon_values, desc='Unpacking epsilon value'):
                epsilon_path = model_path + '/epsilon value = ' + str(epsilon_value)
                deterministic_file = epsilon_path + '/Detereministic data.npy'
                deterministic_data_part = np.load(deterministic_file)
                deteministic_data.append(deterministic_data_part)
                stochastic_file = epsilon_path + '/Stochastic data.npy'
                stochastic_data_part = np.load(stochastic_file)
                stochastic_data.append(stochastic_data_part)

            print('Merging datasets.')
            merged_deteministic_data = np.array(deteministic_data)
            merged_stochastic_data = np.array(stochastic_data)
            print('Saving datasets.')
            np.save(model_path + '/' + 'Detereministic data', merged_deteministic_data)
            np.save(model_path + '/' + 'Stochastic data', merged_stochastic_data)


    #%% Placing merged data into dataframes for plotting figures

    iterations = 100
    N = 1e3
    iterations = 100
    low = 0
    step = 1
    ybins = 50
    stoch_name = 'Stochastic data'
    save_path = os.getcwd() # Path for saving figures and csv files
    normal_person_classes = ['Susceptible',
                             'Singly Infected Sensitive', 'Singly Infected Resistant'
                             ]
    super_inf_person_classes = normal_person_classes + ['Superinfected']
    full_co_inf_person_classes = super_inf_person_classes + ['Double Infected Sensitive', 'Double Infected Resistant']
    default_param_values = {'N': N, 'beta_W': 0.04, 'beta_Z': 0.015, 'gamma': 0.01, 'gamma_T': 0.1}
    model_info = {'Replacement Infection':
                      {'Class Names': normal_person_classes,
                       'Parameter Values': default_param_values},
                  'Exclusive Infection':
                      {'Class Names': normal_person_classes,
                       'Parameter Values': default_param_values},
                  'Unidirectional Conversion':
                      {'Class Names': normal_person_classes,
                       'Parameter Values': {**default_param_values, 'rho': 0.5, 'phi': 0.0}},
                  'Bidirectional Conversion':
                      {'Class Names': normal_person_classes,
                       'Parameter Values': {**default_param_values, 'rho': 0.5, 'phi': 0.05}},
                  'Superinfection':
                      {'Class Names': super_inf_person_classes,
                       'Parameter Values': {**default_param_values, 'q': 0.5}},
                  'Full Coinfection':
                      {'Class Names': full_co_inf_person_classes,
                       'Parameter Values': {**default_param_values, 'q': 0.5}}
                  }

    scenario_info = {'Scenario 1':
                         {'High Epsilon Value': 33},
                     'Scenario 2':
                         {'High Epsilon Value': 50},
                     }

    stochastic_data = {'Scenario': [],
                       'Model': [],
                       'Iteration': [],
                       '% Receiving Treatment': [],
                       'Class': [],
                       'Individuals': []
                       }
    iteration_numbers = [*range(1, iterations+1)]

    for scenario, scenario_info_sub_dict in tqdm.tqdm(scenario_info.items(),desc='Unpacking data for scenario'):
        high = scenario_info_sub_dict['High Epsilon Value']
        epsilon_vals = np.arange(low, high + step, step)
        epsilon_vals_stoch = np.repeat(epsilon_vals, iterations).tolist()
        model_stoch_repeats = len(epsilon_vals_stoch)
        scenario_path = data_directory + '/' + scenario

        for model_name, model_info_sub_dict in tqdm.tqdm(model_info.items(),desc='Unpacking data for model'):
            full_path = scenario_path + '/' + model_name + ' model'
            if not os.path.isfile(save_path + '/stochastic data.csv'):
                npy_name = full_path + '/' + stoch_name + '.npy'
                stoch_data = np.load(npy_name)
                stoch_data_for_time = stoch_data[:, - 1, :, :]
                all_infecteds = {'All Sensitive Infections': [],
                                 'All Resistant Infections': []}
                for index, person_class in enumerate(model_info_sub_dict['Class Names']):
                    value = stoch_data_for_time[:, index, :]
                    value = value.ravel()
                    stochastic_data['Individuals'] += value.tolist()
                    stochastic_data['Class'] += len(value) * [person_class]
                    stochastic_data['% Receiving Treatment'] += epsilon_vals_stoch
                    stochastic_data['Model'] += model_stoch_repeats * [model_name]
                    stochastic_data['Scenario'] += model_stoch_repeats * [scenario]
                    stochastic_data['Iteration'] += len(epsilon_vals) * iteration_numbers
                    if person_class in ['Singly Infected Sensitive', 'Superinfected', 'Double Infected Sensitive']:
                        all_infecteds['All Sensitive Infections'].append(value)
                    if person_class in ['Singly Infected Resistant', 'Superinfected', 'Double Infected Resistant']:
                        all_infecteds['All Resistant Infections'].append(value)
                for person_class, to_sum in all_infecteds.items():
                    value = np.sum(to_sum, axis=0)
                    stochastic_data['Individuals'] += value.tolist()
                    stochastic_data['Class'] += len(value) * [person_class]
                    stochastic_data['% Receiving Treatment'] += epsilon_vals_stoch
                    stochastic_data['Model'] += model_stoch_repeats * [model_name]
                    stochastic_data['Scenario'] += model_stoch_repeats * [scenario]
                    stochastic_data['Iteration'] += len(epsilon_vals) * iteration_numbers

    stochastic_df = pd.DataFrame(stochastic_data)
    stochastic_df.to_csv(save_path+'/stochastic data.csv', index=False)
    # Note you need to run SMC Simulating models deterministically.ipynb in order to have the deterministic data for
    # comparisons. Otherwise, the csv file referenced in the line below will not be there.
    deterministic_df = pd.read_csv(save_path + '/deterministic data.csv')
    deterministic_df['Scenario'] = 'Deterministic Simulation'


    #%% Producing 2D histograms/heatmaps
    plt.figure(figsize=(7.5, 15))
    deterministic_color = 'gold'
    equilibria_color = 'black'
    heatmap_color_map = "coolwarm"
    aspect_ratio = 2
    linewidths = 1.5
    line_plot_steps = 0.25

    model_types = {
        'Conversion models': {
            'model names': ['Unidirectional Conversion', 'Bidirectional Conversion'],
            'model graph titles': ['Unidirectional\nconversion', 'Bidirectional\nconversion']
        },
        'Singly infected models': {
            'model names': ['Exclusive Infection', 'Replacement Infection'],
            'model graph titles': ['Exclusive\nInfection', 'Replacement\nInfection']
        },
        'Double infected models': {
            'model names': ['Superinfection', 'Full Coinfection'],
            'model graph titles': ['Superinfection', 'Full\nco-infection']
        }
    }
    ploting_person_classes = {'names': ['Susceptible', 'All Sensitive Infections', 'All Resistant Infections'],
                              'graph titles': ['Susceptible', 'Infected\nSensitive', 'Infected\nResistant']
                              }

    for scenario, scenario_info_sub_dict in tqdm.tqdm(scenario_info.items(),
                                                      desc='Ploting Graphs for scenario'):
        high = scenario_info_sub_dict['High Epsilon Value']
        xbins = int(high + step)
        bins = [xbins, ybins]
        epsilon_values_as_percent = np.arange(low, high + line_plot_steps, line_plot_steps)
        exclusive_replacement_ees = possible_exclusive_replacement_ee(**default_param_values,
                                                                      epsilon_values_as_percent=epsilon_values_as_percent)
        for model_type, model_type_info in tqdm.tqdm(model_types.items(),
                                                     desc='Ploting Graphs for model types'):
            model_selection = model_type_info['model names']
            model_graph_titles = model_type_info['model graph titles']

            select_deterministic_df = deterministic_df[(deterministic_df.Model.isin(model_selection)) &
                                                       (deterministic_df['% Receiving Treatment'] <= high) &
                                                       (deterministic_df.Class.isin(ploting_person_classes['names']))]
            select_stochastic_df = stochastic_df[(stochastic_df.Model.isin(model_selection)) &
                                                 (stochastic_df.Scenario == scenario) &
                                                 (stochastic_df.Class.isin(ploting_person_classes['names']))]

            fg = sns.FacetGrid(select_stochastic_df,
                               row="Class", col="Model",
                               margin_titles=False, aspect=aspect_ratio)
            # create common colorbar axis
            cax = fg.fig.add_axes([.92, .12, .02, .8])
            vmin = 1
            vmax = iterations
            cbar_ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
            fg.map(sns.histplot, '% Receiving Treatment', 'Individuals', stat="count",
                   bins=bins, vmin=vmin, vmax=vmax,
                   cbar=True, cbar_ax=cax, data=select_stochastic_df, cmap=heatmap_color_map,
                   cbar_kws={"ticks": cbar_ticks})
            for row_index, person_class in enumerate(ploting_person_classes['names']):
                for col_index, model in enumerate(model_selection):
                    ax = fg.axes[row_index, col_index]
                    deterministic_vals = select_deterministic_df[(select_deterministic_df.Model == model) &
                                                                 (select_deterministic_df.Class == person_class)]
                    ax.plot(deterministic_vals['% Receiving Treatment'], deterministic_vals['Individuals'],
                            linewidth=linewidths, color=deterministic_color)
                    if model_type in ['Singly infected models', 'Double infected models']:
                        ee_1, ee_2 = exclusive_replacement_ees
                    else:
                        ee_1, ee_2 = possible_bi_directional_ee(**model_info[model]['Parameter Values'],
                                                                epsilon_values_as_percent=epsilon_values_as_percent)
                    ax.plot(ee_1['% Receiving Treatment'], ee_1[person_class],
                            linewidth=linewidths, color=equilibria_color, linestyle=':')
                    ax.plot(ee_2['% Receiving Treatment'], ee_2[person_class],
                            linewidth=linewidths, color=equilibria_color, linestyle='--')

            fg.set_titles(template='')

            for ax, m in zip(fg.axes[0, :], model_graph_titles):
                ax.set_title(m)
            if model_type == 'Double infected models':
                graph_titles = ['Susceptible', 'All Infected\nWith Sensitive', 'All Infected\nWith Resistant']
            else:
                graph_titles = ploting_person_classes['graph titles']
            for ax, l in zip(fg.axes[:, 0], graph_titles):
                ax.set_ylabel(l)
            fg.set(ylim=(-50, N))
            plt.tight_layout()  # automatically adjust spacing of axes
            fg.fig.subplots_adjust(right=0.9)

            plt.savefig(scenario + ' ' + model_type + '.png')

    #%% Sorting data for production of stacked histograms

    stochastic_pivot_df = stochastic_df.pivot(index=['Scenario', 'Model', 'Iteration', '% Receiving Treatment'],
                                              columns='Class', values='Individuals')
    stochastic_pivot_df.reset_index(inplace=True)

    stochastic_pivot_df['Overall Outcome'] = None
    stochastic_pivot_df.loc[(stochastic_pivot_df['All Sensitive Infections'] == 0) &
                            (stochastic_pivot_df['All Resistant Infections'] == 0), 'Overall Outcome'] = 'Both Extinct'
    stochastic_pivot_df.loc[(stochastic_pivot_df['All Sensitive Infections'] > 0) &
                            (stochastic_pivot_df['All Resistant Infections'] > 0), 'Overall Outcome'] = 'Coexistence'
    stochastic_pivot_df.loc[(stochastic_pivot_df['All Sensitive Infections'] > 0) &
                            (stochastic_pivot_df[
                                 'All Resistant Infections'] == 0), 'Overall Outcome'] = 'Only Sensitive'
    stochastic_pivot_df.loc[(stochastic_pivot_df['All Sensitive Infections'] == 0) &
                            (stochastic_pivot_df['All Resistant Infections'] > 0), 'Overall Outcome'] = 'Only Resistant'

    deterministic_df.Individuals = deterministic_df.Individuals.round()
    deterministic_df.Individuals = deterministic_df.Individuals.clip(lower=0)
    deterministic_pivot_df = deterministic_df.pivot(index=['Scenario', 'Model', '% Receiving Treatment'],
                                                    columns='Class', values='Individuals')
    deterministic_pivot_df.reset_index(inplace=True)
    deterministic_pivot_df = deterministic_pivot_df[
        deterministic_pivot_df['% Receiving Treatment'].isin(np.arange(51, dtype=float).tolist())]

    deterministic_pivot_df['Overall Outcome'] = None
    deterministic_pivot_df.loc[(deterministic_pivot_df['All Sensitive Infections'] == 0) &
                               (deterministic_pivot_df[
                                    'All Resistant Infections'] == 0), 'Overall Outcome'] = 'Both Extinct'
    deterministic_pivot_df.loc[(deterministic_pivot_df['All Sensitive Infections'] > 0) &
                               (deterministic_pivot_df[
                                    'All Resistant Infections'] > 0), 'Overall Outcome'] = 'Coexistence'
    deterministic_pivot_df.loc[(deterministic_pivot_df['All Sensitive Infections'] > 0) &
                               (deterministic_pivot_df[
                                    'All Resistant Infections'] == 0), 'Overall Outcome'] = 'Only Sensitive'
    deterministic_pivot_df.loc[(deterministic_pivot_df['All Sensitive Infections'] == 0) &
                               (deterministic_pivot_df[
                                    'All Resistant Infections'] > 0), 'Overall Outcome'] = 'Only Resistant'
    deterministic_pivot_df['Overall Outcome'].unique()

    pivot_df = pd.concat([stochastic_pivot_df, deterministic_pivot_df], axis=0)

    scenario_info['Deterministic Simulation'] = {'High Epsilon Value': 50}
    #%% Producing stacked histograms

    palette = ['tab:red', 'orange', 'yellow', 'tab:green']
    hue_order = ['Only Resistant', 'Coexistence', 'Only Sensitive', 'Both Extinct']
    for model_type, model_type_info in tqdm.tqdm(model_types.items(),
                                                 desc='Ploting Graphs for model types'):
        selected_data = pivot_df[pivot_df.Model.isin(model_type_info['model names'])]

        fg = sns.FacetGrid(selected_data, row='Scenario', col='Model', sharey=False, aspect=aspect_ratio)
        for row_index, (scenario, scenario_info_sub_dict) in enumerate(scenario_info.items()):
            high = scenario_info_sub_dict['High Epsilon Value']
            bins = int(high + step)
            for col_index, model in enumerate(model_type_info['model names']):
                ax = fg.axes[row_index, col_index]
                selected_data = pivot_df[(pivot_df.Model == model) &
                                         (pivot_df.Scenario == scenario)]
                if row_index == 0 and col_index == 1:
                    legend_inc = True
                else:
                    legend_inc = False
                sns.histplot(selected_data, bins=bins, ax=ax, x='% Receiving Treatment',
                             hue='Overall Outcome', hue_order=hue_order, legend=legend_inc,
                             palette=palette, multiple="stack", linewidth=.5)

        for ax, m in zip(fg.axes[:, 0], ['Scenario 1', 'Scenario 2', 'Deterministic Simulation']):
            ax.set_ylabel(m + '\nCount')
        for row in range(3):
            fg.axes[row, 1].set_ylabel(None)
            fg.axes[row, 1].set_yticklabels([])
        for ax, l in zip(fg.axes[0, :], model_type_info['model graph titles']):
            ax.set_title(l)
        plt.tight_layout()  # automatically adjust spacing of axes

        plt.savefig('Comparing Scenario Outcomes for ' + model_type + '.png')

