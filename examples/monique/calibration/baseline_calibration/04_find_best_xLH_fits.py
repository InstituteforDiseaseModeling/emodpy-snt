# find the xLH value for each admin-seed combination that best matches the reference DHS data
import os
import copy
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from idmtools_calibra.utilities import ll_calculators
from snt.load_paths import load_box_paths

mpl.rcParams['pdf.fonttype'] = 42

# Specify data location

# USER_PATH = r"..\data"
# _, project_path = load_box_paths(user_path=USER_PATH, country_name='Example')
_, project_path = load_box_paths(country_name='Example')


def load_ref_data():
    dhs_pfpr_fname = os.path.join(project_path, 'simulation_inputs', 'prevalence', 'DHS_admin_monthly_microscopy.csv')
    dhs_df = pd.read_csv(dhs_pfpr_fname)
    return dhs_df


def save_best_hab(all_df, working_dir, habitats_fname, create_plot=False):
    all_dhs_df = load_ref_data()
    mult_param = 'Habitat_Multiplier'
    all_df['sim_num_pos'] = all_df['Pop'] * all_df['PfPR U5']
    all_df['sim_pop'] = all_df['Pop'].astype(int)
    all_df['sim_num_pos'] = all_df['sim_num_pos'].astype(int)
    all_admins = sorted(all_df['admin_name'].unique())
    all_runs = sorted(all_df['Run_Number'].unique())
    if create_plot:
        fig1, axes1 = plt.subplots(nrows=len(all_admins), ncols=len(all_runs))
        fig2, axes2 = plt.subplots(nrows=len(all_admins), ncols=len(all_runs))
        fig1.set_figwidth(4 * len(all_runs))
        fig1.set_figheight(4 * len(all_admins))
        fig2.set_figwidth(4 * len(all_runs))
        fig2.set_figheight(4 * len(all_admins))

    hab_df = pd.DataFrame()
    for hfca_run, df in all_df.groupby(['admin_name', 'Run_Number']):
        hfca = hfca_run[0]
        run_number = hfca_run[1]
        row = [i for i, x in enumerate(all_admins) if x == hfca][0]
        col = [i for i, x in enumerate(all_runs) if x == run_number][0]
        dhs_df, df, score_df = return_all_df_vars(mult_param, all_dhs_df, hfca, df, run_number)
        if create_plot:
            # add plot
            fig1, axes1, fig2, axes2, best_LH = plot_graphs(mult_param, score_df, df, dhs_df, fig1, axes1, fig2, axes2,
                                                            row, col, hfca, run_number)
            plt.close()
            fig1.subplots_adjust(wspace=0.5, hspace=0.5)
            fig2.subplots_adjust(wspace=0.5, hspace=0.5)
        else:
            best_LH = get_best_LH(mult_param, score_df)
        hdf = pd.DataFrame({'admin_name': [hfca],
                            'Run_Number': [run_number],
                            'Habitat_Multiplier': [best_LH]})
        hab_df = pd.concat([hab_df, hdf])
    hab_df.to_csv(habitats_fname, index=False)

    if create_plot:
        save_plots(working_dir, fig1, name='loglikelihoods')
        save_plots(working_dir, fig2, name='timeseries')


# save separate plots for each State, instead of all admins in the country in the same file
def plot_outputs_separate(all_df, working_dir):
    all_dhs_df = load_ref_data()
    mult_param = 'Habitat_Multiplier'
    all_df['sim_num_pos'] = all_df['Pop'] * all_df['PfPR U5']
    all_df['sim_pop'] = all_df['Pop'].astype(int)
    all_df['sim_num_pos'] = all_df['sim_num_pos'].astype(int)
    region_names = np.unique(all_dhs_df['NOMREGION'].tolist())
    for rr in range(len(region_names)):
        cur_region = region_names[rr]
        cur_dhs_df = all_dhs_df[all_dhs_df['NOMREGION'] == cur_region]
        all_admins = sorted(np.unique(cur_dhs_df['admin_name'].tolist()))
        sim_df_cur = all_df.loc[all_df['admin_name'].isin(all_admins)]
        all_runs = sorted(sim_df_cur['Run_Number'].unique())
        fig1, axes1 = plt.subplots(nrows=len(all_admins), ncols=len(all_runs))
        fig2, axes2 = plt.subplots(nrows=len(all_admins), ncols=len(all_runs))
        fig1.set_figwidth(4 * len(all_runs))
        fig1.set_figheight(4 * len(all_admins))
        fig2.set_figwidth(4 * len(all_runs))
        fig2.set_figheight(4 * len(all_admins))

        hab_df = pd.DataFrame()
        for hfca_run, df in sim_df_cur.groupby(['admin_name', 'Run_Number']):
            hfca = hfca_run[0]
            run_number = hfca_run[1]
            row = [i for i, x in enumerate(all_admins) if x == hfca][0]
            col = [i for i, x in enumerate(all_runs) if x == run_number][0]
            dhs_df, df, score_df = return_all_df_vars(mult_param=mult_param, all_dhs_df=cur_dhs_df, hfca=hfca, df=df,
                                                      run_number=run_number)
            # add plot
            fig1, axes1, fig2, axes2, best_LH = plot_graphs(mult_param, score_df, df, dhs_df, fig1, axes1, fig2, axes2,
                                                            row, col, hfca, run_number)
            hdf = pd.DataFrame({'admin_name': [hfca],
                                'Run_Number': [run_number],
                                'Habitat_Multiplier': [best_LH]})
            hab_df = pd.concat([hab_df, hdf])

            plt.close()
            fig1.subplots_adjust(wspace=0.5, hspace=0.5)
            fig2.subplots_adjust(wspace=0.5, hspace=0.5)

        save_plots(working_dir, fig1, name='loglikelihoods_%s' % cur_region)
        save_plots(working_dir, fig2, name='timeseries_%s' % cur_region)


# finds dhs_df, df, score_df variables for finalize function
def return_all_df_vars(mult_param, all_dhs_df, hfca, df, run_number):
    dhs_df = copy.copy(all_dhs_df[all_dhs_df['admin_name'] == hfca])
    if len(dhs_df) < 1:
        print(hfca)
        exit()

    scores = []
    for var, sdf in df.groupby(mult_param):
        # merge simulation output into dhs dataframe (for matching month, year values)
        dhs_sim = dhs_df.merge(sdf, on=['admin_name', 'month', 'year'])
        score = np.sum(
            [ll_calculators.beta_binomial(raw_nobs=[x1], sim_nobs=[x2], raw_data=[x3], sim_data=[x4]) for
             x1, x2, x3, x4 in zip(dhs_sim['num_tested'].values,
                                   dhs_sim['sim_pop'].values,
                                   dhs_sim['num_pos'].values,
                                   dhs_sim['sim_num_pos'].values)])

        scores.append(score)
    score_df = pd.DataFrame({mult_param: [var for var, sdf in df.groupby(mult_param)],
                             'score': scores})
    score_df.to_csv(os.path.join(working_dir, 'LL_allxLH_each_admin_run', '%s_run%i.csv' % (hfca, run_number)),
                    index=False)
    return dhs_df, df, score_df


# get the LH that corresponds to the best score
def get_best_LH(mult_param, score_df):
    max_score = np.max(score_df['score'])
    best_LH = score_df[score_df['score'] == max_score][mult_param].values[0]
    return best_LH


# plots all of the graphs for every pair of graphs that represent a region
def plot_graphs(mult_param, score_df, df, dhs_df, fig1, axes1, fig2, axes2, row, col, hfca, run_number):
    df = df.assign(
        date=[datetime.datetime.strptime('15-%i-%i' % (row['month'], row['year']), "%d-%m-%Y").date() for i, row in
              df.iterrows()])
    df = df.sort_values(by=['date', mult_param])
    dhs_df = dhs_df.assign(
        date=[datetime.datetime.strptime('15-%i-%i' % (row['month'], row['year']), "%d-%m-%Y").date() for i, row in
              dhs_df.iterrows()])

    if len(axes1.shape) > 1:
        axes1[row, col].plot(score_df[mult_param], score_df['score'], '-o')
        axes1[row, col].set_xlabel(mult_param)
        axes1[row, col].set_ylabel('log likelihood of DHS U5 PfPR')
        axes1[row, col].set_xscale('log')

        max_score = np.max(score_df['score'])
        best_LH = score_df[score_df['score'] == max_score][mult_param].values[0]
        axes1[row, col].plot(best_LH, max_score, marker='o', markersize=7, color="red")
        axes1[row, col].set_title('%s - run %i\n (best = %.2f)' % (hfca, run_number, best_LH))

        for var, sdf in df.groupby(mult_param):
            if var == best_LH:
                continue
            axes2[row, col].plot(sdf['date'], sdf['PfPR U5'], '-r', linewidth=0.5, alpha=0.3)
        sdf = df[df[mult_param] == best_LH]
        axes2[row, col].plot(sdf['date'], sdf['PfPR U5'], '-r', label='sim')
        axes2[row, col].scatter(dhs_df['date'], [row['num_pos'] / row['num_tested'] for i, row in dhs_df.iterrows()],
                                dhs_df['num_tested'], 'k', label='DHS')
        axes2[row, col].set_xlabel('date')
        axes2[row, col].set_ylabel('U5 PfPR')
        axes2[row, col].set_title('%s - run %i\n (best = %.1f)' % (hfca, run_number, best_LH))
        for label in axes2[row, col].get_xaxis().get_ticklabels()[::2]:
            label.set_visible(False)
    else:
        axes1[row].plot(score_df[mult_param], score_df['score'], '-o')
        axes1[row].set_xlabel(mult_param)
        axes1[row].set_ylabel('log likelihood of DHS U5 PfPR')
        axes1[row].set_xscale('log')

        max_score = np.max(score_df['score'])
        best_LH = score_df[score_df['score'] == max_score][mult_param].values[0]
        axes1[row].plot(best_LH, max_score, marker='o', markersize=7, color="red")
        axes1[row].set_title('%s - run %i\n (best = %.2f)' % (hfca, run_number, best_LH))

        for var, sdf in df.groupby(mult_param):
            if var == best_LH:
                continue
            axes2[row].plot(sdf['date'], sdf['PfPR U5'], '-r', linewidth=0.5, alpha=0.3)
        sdf = df[df[mult_param] == best_LH]
        axes2[row].plot(sdf['date'], sdf['PfPR U5'], '-r', label='sim')
        axes2[row].scatter(dhs_df['date'], [row['num_pos'] / row['num_tested'] for i, row in dhs_df.iterrows()],
                           dhs_df['num_tested'], 'k', label='DHS')
        axes2[row].set_xlabel('date')
        axes2[row].set_ylabel('U5 PfPR')
        axes2[row].set_title('%s - run %i\n (best = %.1f)' % (hfca, run_number, best_LH))
        for label in axes2[row].get_xaxis().get_ticklabels()[::2]:
            label.set_visible(False)
    return fig1, axes1, fig2, axes2, best_LH


# saves plots
def save_plots(working_dir, fig, name):
    fig.savefig(os.path.join(working_dir, 'compare_xLH_sweep_with_DHS_%s.png' % name), bbox_inches='tight',
                format='PNG')


if __name__ == "__main__":

    expt_name = 'baseline_calibration'
    working_dir = os.path.join(project_path, 'simulation_output', 'baseline_calibration', expt_name)
    if not os.path.exists(os.path.join(working_dir, 'LL_allxLH_each_admin_run')):
        os.mkdir(os.path.join(working_dir, 'LL_allxLH_each_admin_run'))

    habitats_fname = os.path.join(project_path, 'simulation_inputs', 'larval_habitats',
                                  'larval_habitat_multipliers_v1.csv')

    all_df = pd.read_csv(os.path.join(working_dir, 'monthly_U5_PfPR.csv'))

    save_best_hab(all_df, working_dir, habitats_fname, create_plot=False)
    plot_outputs_separate(all_df, working_dir)
