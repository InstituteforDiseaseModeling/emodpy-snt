import os
import datetime
import pandas as pd
import numpy as np
from idmtools.entities import IAnalyzer


class monthlyEventAnalyzerITN(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026, output_file_suffix=''):
        super(monthlyEventAnalyzerITN, self).__init__(working_dir=working_dir,
                                                   filenames=["output/ReportEventCounter.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment',  # 'Received_Self_Medication',
                             'Bednet_Using', # 'Bednet_Got_New_One',
                             # currently removed 'Bednet_Got_New_One', since length is 1 longer than expected for unknown reasons
                             'Received_Campaign_Drugs', 'Received_IRS', 'Received_Vaccine', 'Received_PMC_VaccDrug']
        else:
            self.channels = channels
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.output_file_suffix = output_file_suffix

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        channels_in_expt = [x for x in self.channels if x in data[self.filenames[0]]['Channels'].keys()]

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in channels_in_expt})
        simdata['Time'] = simdata.index + 1

        simdata['Day'] = simdata['Time'] % 365
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)

        for missing_channel in [x for x in self.channels if x not in channels_in_expt]:
            simdata[missing_channel] = 0

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True, )
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        df = adf.groupby(['admin_name', 'date', 'Run_Number'])[self.channels].agg(np.sum).reset_index()
        df.to_csv(os.path.join(self.working_dir, self.expt_name, 'monthly_Event_Count%s.csv' % self.output_file_suffix), index=False)



class monthlyEventAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026, output_file_suffix=''):
        super(monthlyEventAnalyzer, self).__init__(working_dir=working_dir,
                                                   filenames=["output/ReportEventCounter.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment',
                             'Bednet_Using', # 'Bednet_Got_New_One',  # 'Received_Self_Medication',
                             # currently removed 'Bednet_Got_New_One', since length is 1 longer than expected for unknown reasons
                             'Received_Campaign_Drugs', 'Received_IRS', 'Received_Vaccine', 'Received_PMC_VaccDrug']
        else:
            self.channels = channels
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.output_file_suffix = output_file_suffix

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        channels_in_expt = [x for x in self.channels if x in data[self.filenames[0]]['Channels'].keys()]

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in channels_in_expt})
        simdata['Time'] = simdata.index

        simdata['Day'] = simdata['Time'] % 365
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)

        for missing_channel in [x for x in self.channels if x not in channels_in_expt]:
            simdata[missing_channel] = 0

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True, )
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        df = adf.groupby(['admin_name', 'date', 'Run_Number'])[self.channels].agg(np.sum).reset_index()
        df.to_csv(os.path.join(self.working_dir, self.expt_name, 'monthly_Event_Count%s.csv' % self.output_file_suffix), index=False)


class monthlyUsageLLIN(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyUsageLLIN, self).__init__(working_dir=working_dir,
                                               filenames=["output/ReportEventCounter.json",
                                                          "output/ReportMalariaFiltered.json"]
                                               )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Bednet_Using']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # added to bypass failed cases
    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        channels_in_expt = [x for x in self.channels if x in data[self.filenames[0]]['Channels'].keys()]
        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in channels_in_expt})
        simdata['Time'] = simdata.index

        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
        d['Time'] = d.index

        if len(channels_in_expt) > 0:
            simdata = pd.merge(left=simdata, right=d, on='Time')
        else:
            simdata = d
        for missing_channel in [x for x in self.channels if x not in channels_in_expt]:
            simdata[missing_channel] = 0

        simdata['day_of_year'] = simdata['Time'] % 365
        simdata['month'] = simdata['day_of_year'].apply(lambda x: self.monthparser((x + 1) % 365))
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        mean_channels = self.channels + ['Statistical Population']
        adf = adf.groupby(['admin_name', 'date', 'Run_Number'])[mean_channels].agg(np.mean).reset_index()
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'MonthlyUsageLLIN.csv'), index=False)



if __name__ == "__main__":
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform
    from snt.load_paths import load_box_paths

    platform = Platform('Calculon')

    data_path, project_path = load_box_paths(country_name='Burundi')

    working_dir = os.path.join(project_path, 'simulation_output', '2010_to_present')
    start_year = 2010  # simulation starts in January of this year
    end_year = 2021  # simulation ends in December of this year
    # start_year = 2021  # simulation starts in January of this year
    # end_year = 2030  # simulation ends in December of this year

    expt_ids = {
        'test_analyzers_from_toPresent_v3': '5c0726d3-4cf3-ed11-aa06-b88303911bc1'
        # 'NGA_toPresent_allInter': '0a297502-088c-ed11-aa00-b88303911bc1',
    }
    include_LLINp = False  # determines whether number of new infections among individuals with/without LLINps obtained
    itn_comparison = False

    if (not include_LLINp) and (not itn_comparison):
        for expname, expid in expt_ids.items():
            print('running expt %s' % expname)
            report_count_channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment',
                                     'Received_Self_Medication', 'Bednet_Got_New_One', 'Bednet_Using',
                                     'Received_Campaign_Drugs', 'Received_IRS'
                                     ]
            report_count_channels = None



            analyzers = [

                monthlyEventAnalyzer(expt_name=expname,
                                     channels=report_count_channels,
                                     sweep_variables=["Run_Number", "admin_name"],
                                     working_dir=working_dir,
                                     start_year=start_year,
                                     end_year=end_year),


            ]
            am = AnalyzeManager(platform=platform, ids=[(expid, ItemType.EXPERIMENT)], analyzers=analyzers,
                                analyze_failed_items=True)
            am.analyze()

    elif include_LLINp:
        for expname, expid in expt_ids.items():
            print('running expt %s' % expname)
            analyzers = [

                monthlyEventAnalyzer(expt_name=expname,
                                     sweep_variables=["Run_Number", "admin_name"],
                                     working_dir=working_dir,
                                     start_year=start_year,
                                     end_year=end_year),

            ]
            am = AnalyzeManager(platform=platform, ids=[expid, ItemType.EXPERIMENT], analyzers=analyzers,
                                force_analyze=True)
            am.analyze()

