import os
import datetime
import logging
import pandas as pd
import numpy as np
from idmtools.entities import IAnalyzer

logger = logging.getLogger(__name__)


class seasonalityIncidenceAnalyzer(IAnalyzer):
    """
    Get incidence across months from simulations (to be compared against reference datasets)
    """

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, working_dir="."):
        super().__init__(working_dir=working_dir, filenames=['output/ReportEventCounter.json',
                                                             'output/ReportMalariaFiltered.json'])

        self.population_channel = 'Statistical Population'
        self.case_channel = 'Received_Treatment'
        self.prev_channel = 'PfHRP2 Prevalence'
        self.nmf_channel = 'Received_NMF_Treatment'
        self.comparison_channel = 'Treated Cases NMF Adjusted'
        self.working_dir = working_dir

    def map(self, data, simulation):
        """
        Extract data from output data and accumulate in same bins as reference.
        """

        # Load data from simulation
        simdata = {self.case_channel: data[self.filenames[0]]['Channels'][self.case_channel]['Data'][-365:],
                   self.nmf_channel: data[self.filenames[0]]['Channels'][self.nmf_channel]['Data'][-365:],
                   self.population_channel: data[self.filenames[1]]['Channels'][self.population_channel]['Data'][-365:],
                   self.prev_channel: data[self.filenames[1]]['Channels'][self.prev_channel]['Data'][-365:]}

        simdata = pd.DataFrame(simdata)
        simdata[self.comparison_channel] = simdata[self.case_channel] + simdata[self.nmf_channel]

        simdata = simdata[-365:].reset_index(drop=True)
        simdata['Time'] = simdata.index
        simdata['Day'] = simdata['Time'] % 365
        simdata['Month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))

        simdata = simdata.rename(columns={self.population_channel: 'Trials',
                                          self.comparison_channel: 'Observations'})

        s1 = simdata.groupby('Month')['Trials'].agg(np.mean).reset_index()
        s2 = simdata.groupby('Month')['Observations'].agg(np.sum).reset_index()
        simdata = pd.merge(left=s1, right=s2, on='Month')
        simdata = simdata[['Month', 'Trials', 'Observations']]
        simdata['CasesPer1000'] = simdata['Observations'] / simdata['Trials'] * 1000

        # add in the simulation id for debugging
        # simdata['sim_id'] = simulation.id.hex

        # simdata = simdata.set_index(['Month'])

        return simdata

    def reduce(self, all_data):
        """
        Calculate the mean output result for each experiment. Note that we assume the setup is that each experiment
            uses a single set of parameters and we compare its mean incidence against the reference.
        """
        # selected = list(all_data.values())
        #
        # # Stack selected_data from each parser, adding unique (sim_id) and shared (sample) levels to MultiIndex
        # combine_levels = ['sample', 'sim_id', 'Counts']
        # combined = pd.concat(selected, axis=1,
        #                      keys=[(s.tags.get('__sample_index__'), s.id) for s in all_data.keys()],
        #                      names=combine_levels)
        #
        # data = combined.groupby(level=['sample', 'Counts'], axis=1).mean()
        # compare_results = data.groupby(level='sample', axis=1).apply(self.compare)
        # # Make sure index is sorted in correct order
        # compare_results.index = compare_results.index.astype(int)
        # compare_results = compare_results.sort_index(ascending=True)

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return
        # if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
        #     os.mkdir(os.path.join(self.working_dir, self.expt_name))
        adf = pd.concat(selected).reset_index(drop=True)
        adf = adf.groupby(['Month']).agg(np.mean).reset_index()
        adf.to_csv(os.path.join(self.working_dir, 'All_Age_monthly_prevalence.csv'), index=False)
        return adf





