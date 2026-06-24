import os
import pandas as pd
import numpy as np
import datetime
from idmtools.entities import IAnalyzer


class MonthlyPfPRU5Analyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2010, end_year=2016,
                 input_filename_base='MalariaSummaryReport_Monthly'):
        super(MonthlyPfPRU5Analyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/%s%d.json" % (input_filename_base, x)
                                                               for x in range(start_year, (1 + end_year))]
                                                    )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.mult_param = 'Habitat_Multiplier'
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    def map(self, data, simulation):

        adf = pd.DataFrame()
        aa = 1  # index for the 0.5-5 year age group
        for year, fname in zip(range(self.start_year, (self.end_year + 1)), self.filenames):
            # population size
            pop = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][
                  :12]  # remove final five days: assume final five days have same average as rest of month
            pop_monthly = [x[aa] for x in pop]
            # PfPR
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            pfpr_monthly = [x[aa] for x in d]
            # combine in data frame
            simdata = pd.DataFrame({'month': list(range(1, 13)),
                                    'Pop': pop_monthly,
                                    'PfPR U5': pfpr_monthly,
                                    })
            simdata['year'] = year
            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        all_df = pd.concat(selected).reset_index(drop=True)

        # don't take average across runs; calculate the best xLH for each seed (each of which uses different quantile)
        # all_df = all_df.groupby([self.mult_param, 'month', 'admin_name'])[['PfPR U5', 'Pop']].agg(np.mean).reset_index()
        # all_df = all_df.sort_values(by=[self.mult_param, 'month', 'archetype', 'DS_Name_for_ITN'])
        all_df.to_csv(os.path.join(self.working_dir, 'monthly_U5_PfPR.csv'), index=False)




class monthlyU5PrevalenceAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyU5PrevalenceAnalyzer, self).__init__(working_dir=working_dir,
                                                        filenames=["output/ReportMalariaFiltered__RDT_mic_PfPR_U5.json"]
                                                        )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'True Prevalence',
                               'PCR Parasite Prevalence', 'Blood Smear Parasite Prevalence', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # added to bypass failed cases
    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):
        d = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.inset_channels})
        d['Time'] = d.index
        simdata = d
        simdata['Day'] = simdata['Time'] % 365
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        # adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        sum_channels = ['New Clinical Cases']
        mean_channels = ['Statistical Population', 'True Prevalence', 'PCR Parasite Prevalence',
                         'Blood Smear Parasite Prevalence', 'PfHRP2 Prevalence']

        df = adf.groupby(['year', 'month']+self.sweep_variables)[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(['year', 'month']+self.sweep_variables)[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=['year', 'month']+self.sweep_variables)
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'monthly_U5_PfPR_mic_RDT.csv'), index=False)








if __name__ == "__main__":
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform
    from snt.load_paths import load_box_paths

    platform = Platform('Calculon')

    data_path, project_path = load_box_paths(country_name='Nigeria')

    working_dir = os.path.join(project_path, 'simulation_outputs', 'baseline_calibration')
    start_year = 2010  # simulation starts in January of this year
    end_year = 2021  # simulation ends in December of this year


    expt_ids = {
        'PfPR_sweep_main_NGA_v1_testRDT': '58ddf0c3-90a9-ee11-9eff-b88303912b51'
    }

    for expname, expid in expt_ids.items():
        print('running expt %s' % expname)

        sweep_variables = ["Run_Number",
                           "Habitat_Multiplier",
                           "admin_name",
                           ]
        analyzers = [MonthlyPfPRU5Analyzer(expname, sweep_variables, working_dir, start_year, end_year),
                     monthlyU5PrevalenceAnalyzer(expname, sweep_variables, working_dir, start_year, end_year)]


        am = AnalyzeManager(platform=platform, ids=[(expid, ItemType.EXPERIMENT)], analyzers=analyzers,
                            analyze_failed_items=True)

        am.analyze()
