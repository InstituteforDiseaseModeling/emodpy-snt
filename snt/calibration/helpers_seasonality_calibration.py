import os
import pandas as pd
from snt.helpers_sim_setup import habitat_scales


# fit values for all months
def get_burnin_spline_values(hfca, month_scalar, project_path):
    month_vals = [month_scalar] * 12

    hdf = habitat_scales(project_path)
    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return month_vals, fraction


def get_spline_values(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['Month%d' % x for x in range(1, 13)],
                       'Guess': [0.005] * 12,
                       'Min': [0.00001] * 12,
                       'Max': [0.01] * 12})
    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [11], 'Min': [9], 'Max': [12.5]}),
                    pd.DataFrame({'Name': ['Constant'], 'Guess': [6], 'Min': [6], 'Max': [10]})])

    # df.loc[df['Name'] == 'Month1', 'Guess'] = 0
    # df.loc[df['Name'] == 'Month2', 'Guess'] = 0
    # df.loc[df['Name'] == 'Month3', 'Guess'] = 0
    # df.loc[df['Name'] == 'Month4', 'Guess'] = 0
    # df.loc[df['Name'] == 'Month5', 'Guess'] = 0.001
    # df.loc[df['Name'] == 'Month6', 'Guess'] = 0.001
    # df.loc[df['Name'] == 'Month7', 'Guess'] = 0.00775
    # df.loc[df['Name'] == 'Month8', 'Guess'] = 0.0081
    # df.loc[df['Name'] == 'Month9', 'Guess'] = 0.003
    # df.loc[df['Name'] == 'Month10', 'Guess'] = 0.002
    # df.loc[df['Name'] == 'Month11', 'Guess'] = 0
    # df.loc[df['Name'] == 'Month12', 'Guess'] = 0
    #
    # df['Min'] = df['Guess'] * 0.1
    # df['Max'] = df['Guess'] * 10
    df['Dynamic'] = True

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# select num_points months to fit (one of them is always January): will be num_points*2-1 values (num_points - 1) months selected plus num_points values for each month
def get_spline_values2(hfca, project_path, num_points=4):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthNum%d' % x for x in range(1, (num_points + 1))] + ['MonthVal%d' % x for x in
                                                                                         range(1, (num_points + 1))],
                       'Guess': [0] * (num_points * 2),
                       'Min': [0] * (num_points * 2),
                       'Max': [0] * (num_points * 2)})

    for month_num in range(1, (num_points + 1)):
        num_name = 'MonthNum%d' % month_num
        val_name = 'MonthVal%d' % month_num
        df.loc[df['Name'] == num_name, 'Guess'] = month_num
        df.loc[df['Name'] == num_name, 'Min'] = 1.51
        df.loc[df['Name'] == num_name, 'Max'] = 12.49
        df.loc[df['Name'] == val_name, 'Guess'] = 0.005
        df.loc[df['Name'] == val_name, 'Min'] = 0.00001
        df.loc[df['Name'] == val_name, 'Max'] = 0.01  # 10

    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [11], 'Min': [9], 'Max': [12.5]})])
    df['Dynamic'] = True

    # MonthNum1 set to January (necessary for how I've set up the periodic spline)
    df.loc[df['Name'] == 'MonthNum1', 'Guess'] = 1  # January is always set as the first month
    df.loc[df['Name'] == 'MonthNum1', 'Min'] = 0
    df.loc[df['Name'] == 'MonthNum1', 'Dynamic'] = False  # month 1 not fitted

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# fit values for odd months, cubic spline to get values for even months
def get_spline_values3(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthVal%d' % (x * 2 + 1) for x in range(6)],
                       'Guess': [0.0075] * 6,
                       'Min': [0.00001] * 6,
                       'Max': [0.01] * 6})

    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [11], 'Min': [9], 'Max': [12.5]})])
    df['Dynamic'] = True

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# fit values for all months
def get_spline_values4(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthVal%d' % x for x in range(1, 13)],
                       'Guess': [0.01] * 12,
                       'Min': [0.00001] * 12,
                       'Max': [1] * 12})  # !0.01
    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [11], 'Min': [8], 'Max': [12.2]})])
    df['Dynamic'] = True

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# fit values for all months
def get_spline_values4_constantMaxHab(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthVal%d' % x for x in range(1, 13)],
                       'Guess': [0.01] * 12,
                       'Min': [0.00001] * 12,
                       'Max': [1] * 12})  # !0.01
    df['Dynamic'] = True

    df = pd.concat(
        [df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [12], 'Min': [12], 'Max': [12], 'Dynamic': [False]})])

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction



# fit values for max hab scalar, keeping monthly values constant
def get_spline_values4_constantMonthly(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthVal%d' % x for x in range(1, 13)],
                       'Guess': [0.01] * 12,
                       'Min': [0.01] * 12,
                       'Max': [0.01] * 12})
    df['Dynamic'] = False
    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [10], 'Min': [8], 'Max': [12.2], 'Dynamic': [True]})])

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# fit values for all months
def get_spline_values5(hfca, project_path):
    hdf = habitat_scales(project_path)

    df = pd.DataFrame({'Name': ['MonthVal%d' % x for x in range(1, 13)],
                       'Guess': [0.0075] * 12,
                       'Min': [0.00001] * 12,
                       'Max': [0.02] * 12})  # !0.01

    df = pd.concat([df, pd.DataFrame({'Name': ['MaxHab'], 'Guess': [11], 'Min': [9], 'Max': [12.5]})])
    df['Dynamic'] = True

    # for starting guess on monthly multipliers, use rainfall proportional to maximum rainfall month
    rain_fname = os.path.join(project_path, 'SpatialClustering', 'input_layers', 'rain', 'mean_rain_values_by_DS.csv')
    rain_df = pd.read_csv(rain_fname)
    rain_df = rain_df[rain_df['DS'] == hfca].reset_index()
    max_rain = max(
        [rain_df.iloc[0, col] for col in range(len(rain_df.columns)) if 'rainfall_month' in rain_df.columns[col]])
    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        if month_num < 10:
            col_name = 'rainfall_month_0%d' % month_num
        else:
            col_name = 'rainfall_month_%d' % month_num
        rain_col = [col for col in rain_df.columns if col_name in col]
        df.loc[df['Name'] == val_name, 'Guess'] = rain_df.loc[0, rain_col[0]] / max_rain * 0.01

    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0, a), max(0, f), max(0, g))

    return df, fraction


# use the parameter values with the highest likelihood from a previous iteration as the starting point
def update_starting_spline_values(spline, round1_best_df):
    spline.loc[spline['Name'] == 'MaxHab', 'Guess'] = round1_best_df.MaxHab[0]
    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        spline.loc[spline['Name'] == val_name, 'Guess'] = round1_best_df.loc[0, val_name]

    return spline


# use the parameter values with the highest likelihood from a previous iteration as the starting point
def update_starting_spline_values_constantMaxHab(spline, round1_best_df):
    spline.loc[spline['Name'] == 'MaxHab', 'Guess'] = round1_best_df.MaxHab[0]
    spline.loc[spline['Name'] == 'MaxHab', 'Min'] = round1_best_df.MaxHab[0] * 0.98
    spline.loc[spline['Name'] == 'MaxHab', 'Max'] = round1_best_df.MaxHab[0] * 1.02

    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        spline.loc[spline['Name'] == val_name, 'Guess'] = round1_best_df.loc[0, val_name]

    return spline


# use the parameter values with the highest likelihood from a previous iteration as the starting point
def update_starting_spline_values_constantMonthly(spline, round1_best_df):
    spline.loc[spline['Name'] == 'MaxHab', 'Guess'] = round1_best_df.MaxHab[0]

    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        spline.loc[spline['Name'] == val_name, 'Guess'] = round1_best_df.loc[0, val_name]
        spline.loc[spline['Name'] == val_name, 'Min'] = round1_best_df.loc[0, val_name] * 0.98
        spline.loc[spline['Name'] == val_name, 'Max'] = round1_best_df.loc[0, val_name] * 1.02

    return spline


def set_starting_spline_values_from_file(spline, initial_month_values, rep_admin):
    cur_month_values = initial_month_values[initial_month_values['archetype'] == rep_admin]
    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        spline.loc[spline['Name'] == val_name, 'Guess'] = cur_month_values.loc[0, val_name]

    return spline



def set_starting_spline_values_from_file_constantMonthly(spline, initial_month_values, rep_admin):
    cur_month_values = initial_month_values[initial_month_values['archetype'] == rep_admin]
    for month_num in range(1, 13):
        val_name = 'MonthVal%d' % month_num
        spline.loc[spline['Name'] == val_name, 'Guess'] = cur_month_values.iloc[0][val_name]
        spline.loc[spline['Name'] == val_name, 'Min'] = cur_month_values.iloc[0][val_name] * 0.98
        spline.loc[spline['Name'] == val_name, 'Max'] = cur_month_values.iloc[0][val_name] * 1.02

    spline.loc[spline['Name'] == 'MaxHab', 'Guess'] = cur_month_values.iloc[0]['MaxHab']

    return spline

def get_cases(hfca, project_path):
    reference_fname = os.path.join(project_path, 'simulation_inputs', 'incidence', 'archetype_incidence.csv')

    ref_df = pd.read_csv(reference_fname)
    ref_df = ref_df[ref_df['seasonality_archetype'] == hfca]
    ref_df = ref_df.rename(columns={'month': 'Month',
                                    'population': 'Trials'})
    ref_df['Observations'] = ref_df['incidence'] * ref_df['Trials'] / 1000

    return ref_df[['Month', 'Trials', 'Observations']]
