import os
import datetime
import pandas as pd
import numpy as np
from idmtools.entities import IAnalyzer


class monthlyU1PfPRAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyU1PfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/MalariaSummaryReport_Monthly_%d.json" % x
                                                               for x in range(start_year, end_year)]
                                                    )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # def filter(self, simulation):
    #     return simulation.tags["admin_name"] == ''

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
            pfpr = [x[0] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:12]
            clinical_cases = [x[0] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:12]
            severe_cases = [x[0] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]  # this add pop col in U1
            pop = [x[0] for x in d]
            simdata = pd.DataFrame({'month': range(1, 13),
                                    'PfPR U1': pfpr,
                                    'Cases U1': clinical_cases,
                                    'Severe cases U1': severe_cases,
                                    'Pop U1': pop})
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'U1_PfPR_ClinicalIncidence.csv')), index=False)


class monthlyU5PfPRAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyU5PfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/MalariaSummaryReport_Monthly%d.json" % x
                                                               for x in range(start_year, end_year + 1)]
                                                    )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year + 1), self.filenames):
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
            pfpr = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:12]
            clinical_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:12]
            severe_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]
            pop = [x[1] for x in d]
            simdata = pd.DataFrame({'month': range(1, 13),
                                    'PfPR U5': pfpr,
                                    'Cases U5': clinical_cases,
                                    'Severe cases U5': severe_cases,
                                    'Pop U5': pop})
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'U5_PfPR_ClinicalIncidence.csv')), index=False)


class MonthlyPfPRAnalyzerByAge(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(MonthlyPfPRAnalyzerByAge, self).__init__(working_dir=working_dir,
                                                       filenames=["output/MalariaSummaryReport_Monthly_U1U5_%d.json" % x
                                                                  for x in range(start_year, end_year)]  # ,2020
                                                       )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        #        self.agebins = agebins or [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 15, 30, 50, 125]
        self.agebins = [1, 5, 120]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # def filter(self, simulation):
    #     return simulation.tags["admin_name"] == ''

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):  # , 2020
            for age in list(range(0, len(self.agebins))):
                d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
                pfpr = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:12]
                clinical_cases = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:12]
                severe_cases = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]  # this add pop col in U5
                pop = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Mild Anemia by Age Bin'][:12]
                mild_anemia = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Moderate Anemia by Age Bin'][:12]
                moderate_anemia = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Anemia by Age Bin'][:12]
                severe_anemia = [x[age] for x in d]
                simdata = pd.DataFrame({'month': range(1, 13),
                                        'PfPR': pfpr,
                                        'Cases': clinical_cases,
                                        'Severe cases': severe_cases,
                                        'Pop': pop,
                                        'Mild anaemia': mild_anemia,
                                        'Moderate anaemia': moderate_anemia,
                                        'Severe anaemia': severe_anemia
                                        })
                simdata['year'] = year
                simdata['agebin'] = self.agebins[age]
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'Agebins_PfPR_ClinicalIncidenceAnemia.csv')),
                   index=False)


class monthlyTreatedCasesAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                          filenames=["output/ReportEventCounter.json",
                                                                     "output/ReportMalariaFiltered.json"]
                                                          )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment', 'Received_Vaccine']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
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
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        sum_channels = self.channels + ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = adf.groupby(['admin_name', 'date', 'Run_Number'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(['admin_name', 'date', 'Run_Number'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=['admin_name', 'date', 'Run_Number'])
        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_monthly_Cases.csv'), index=False)


class monthlyPrevalenceAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyPrevalenceAnalyzer, self).__init__(working_dir=working_dir,
                                                        filenames=["output/ReportMalariaFiltered.json"]
                                                        )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'True Parasite Prevalence',
                               'PCR Parasite Prevalence', 'Blood Smear Parasite Prevalence', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # added to bypass failed cases
    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):
        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        sum_channels = ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'True Parasite Prevalence', 'PCR Parasite Prevalence',
                         'Blood Smear Parasite Prevalence', 'PfHRP2 Prevalence']

        df = adf.groupby(['admin_name', 'date', 'Run_Number'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(['admin_name', 'date', 'Run_Number'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=['admin_name', 'date', 'Run_Number'])
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_monthly_prevalence.csv'), index=False)


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


class monthlySevereTreatedByAgeAnalyzer(IAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, event_name='Received_Severe_Treatment', agebins=None,
                 sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):

        super(monthlySevereTreatedByAgeAnalyzer, self).__init__(working_dir=working_dir,
                                                                filenames=["output/ReportEventRecorder.csv"]
                                                                )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        self.event_name = event_name
        self.agebins = agebins or [1, 5, 200]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # def filter(self, simulation):
    #     return simulation.tags['admin_name'] == 'Batie'

    def map(self, data, simulation):

        output_data = data[self.filenames[0]]
        output_data = output_data[output_data['Event_Name'] == self.event_name]

        simdata = pd.DataFrame()
        if len(output_data) > 0:  # there are events of this type
            # output_data['Day'] = output_data['Time'] % 365
            # output_data['month'] = output_data['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
            # output_data['year'] = output_data['Time'].apply(lambda x: int(x / 365) + self.start_year)
            # output_data['age in years'] = output_data['Age'] / 365
            output_data = output_data.copy()  # ensure it's not a view

            output_data.loc[:, 'Day'] = output_data['Time'] % 365
            output_data.loc[:, 'month'] = output_data['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
            output_data.loc[:, 'year'] = output_data['Time'] // 365 + self.start_year
            output_data.loc[:, 'age in years'] = output_data['Age'] / 365
            for agemax in self.agebins:
                if agemax < 200:
                    agelabel = 'U%d' % agemax
                else:
                    agelabel = 'all_ages'
                if agemax == 5:
                    agemin = 0.25
                else:
                    agemin = 0
                d = output_data[(output_data['age in years'] < agemax) & (output_data['age in years'] > agemin)]
                g = d.groupby(['year', 'month'])['Event_Name'].agg(len).reset_index()
                g = g.rename(columns={'Event_Name': 'Num_%s_Received_Severe_Treatment' % agelabel})
                if simdata.empty:
                    simdata = g
                else:
                    if not g.empty:
                        simdata = pd.merge(left=simdata, right=g, on=['year', 'month'], how='outer')
                        simdata = simdata.fillna(0)

            for sweep_var in self.sweep_variables:
                if sweep_var in simulation.tags.keys():
                    simdata[sweep_var] = simulation.tags[sweep_var]
        else:
            simdata = pd.DataFrame(columns=['year', 'month', 'Num_U5_Received_Severe_Treatment',
                                            'Num_U1_Received_Severe_Treatment',
                                            'Num_all_ages_Received_Severe_Treatment'] + self.sweep_variables)
        return simdata

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf = adf.fillna(0)
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'Treated_Severe_Monthly_Cases_By_Age.csv'),
                   index=False)

        included_child_bins = ['U%i' % x for x in self.agebins if x < 20]
        for agelabel in included_child_bins:
            severe_treat_df = adf[
                ['year', 'month', 'Num_%s_Received_Severe_Treatment' % agelabel] + self.sweep_variables]
            # cast to int65 data type for merge with incidence df
            severe_treat_df = severe_treat_df.astype({'month': 'int64', 'year': 'int64', 'Run_Number': 'int64'})

            # combine with existing columns of the U5 clinical incidence and PfPR dataframe
            incidence_df = pd.read_csv(
                os.path.join(self.working_dir, self.expt_name, '%s_PfPR_ClinicalIncidence.csv' % agelabel))
            merged_df = pd.merge(left=incidence_df, right=severe_treat_df,
                                 on=['admin_name', 'year', 'month', 'Run_Number'],
                                 how='left')
            merged_df = merged_df.fillna(0)

            # fix any excess treated cases!
            merged_df['num severe cases %s' % agelabel] = merged_df['Severe cases %s' % agelabel] * merged_df[
                'Pop %s' % agelabel] * 30 / 365
            merged_df['excess sev treat %s' % agelabel] = merged_df['Num_%s_Received_Severe_Treatment' % agelabel] - \
                                                          merged_df['num severe cases %s' % agelabel]

            for (rn, admin_name), rdf in merged_df.groupby(['Run_Number', 'admin_name']):
                for r, row in rdf.iterrows():
                    if row['excess sev treat %s' % agelabel] < 1:
                        continue
                    # fix Jan 2020 (start of sim) excess treated severe cases
                    if row['year'] == self.start_year and row['month'] == 1:
                        merged_df.loc[(merged_df['year'] == self.start_year) & (merged_df['month'] == 1) & (
                                merged_df['Run_Number'] == rn) & (merged_df['admin_name'] == admin_name),
                                      'Num_%s_Received_Severe_Treatment' % agelabel] = np.sum(
                            merged_df[(merged_df['year'] == self.start_year) &
                                      (merged_df['month'] == 1) &
                                      (merged_df['Run_Number'] == rn) &
                                      (merged_df['admin_name'] == admin_name)]['num severe cases %s' % agelabel])
                    else:
                        # figure out which is previous month
                        newyear = row['year']
                        newmonth = row['month'] - 1
                        if newmonth < 1:
                            newyear -= 1
                        excess = row['excess sev treat %s' % agelabel]
                        merged_df.loc[(merged_df['year'] == self.start_year) & (merged_df['month'] == 1) & (
                                merged_df['Run_Number'] == rn) & (merged_df[
                                                                      'admin_name'] == admin_name), 'Num_%s_Received_Severe_Treatment' % agelabel] = \
                            merged_df.loc[(merged_df['year'] == self.start_year) & (merged_df['month'] == 1) & (
                                    merged_df['Run_Number'] == rn) & (merged_df['admin_name'] == admin_name),
                                          'Num_%s_Received_Severe_Treatment' % agelabel] - excess
                        merged_df.loc[(merged_df['year'] == self.start_year) & (merged_df['month'] == 1) & (
                                merged_df['Run_Number'] == rn) & (merged_df[
                                                                      'admin_name'] == admin_name), 'Num_%s_Received_Severe_Treatment' % agelabel] = \
                            merged_df.loc[(merged_df['year'] == self.start_year) & (merged_df['month'] == 1) & (
                                    merged_df['Run_Number'] == rn) & (merged_df['admin_name'] == admin_name),
                                          'Num_%s_Received_Severe_Treatment' % agelabel] + excess
            merged_df['excess sev treat %s' % agelabel] = merged_df['Num_%s_Received_Severe_Treatment' % agelabel] - \
                                                          merged_df['num severe cases %s' % agelabel]
            merged_df.loc[
                merged_df['excess sev treat %s' % agelabel] > 0.5, 'Num_%s_Received_Severe_Treatment' % agelabel] = \
                merged_df.loc[merged_df['excess sev treat %s' % agelabel] > 0.5, 'num severe cases %s' % agelabel]

            del merged_df['num severe cases %s' % agelabel]
            del merged_df['excess sev treat %s' % agelabel]
            merged_df.to_csv(os.path.join(self.working_dir, self.expt_name,
                                          '%s_PfPR_ClinicalIncidence_severeTreatment.csv' % agelabel), index=False)


class MonthlyNewInfectionsAnalyzer(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026,
                 input_filename_base='MalariaSummaryReport_Monthly',
                 output_filename='newInfections_PfPR_cases_monthly_byAgeGroup.csv'):

        super(MonthlyNewInfectionsAnalyzer, self).__init__(working_dir=working_dir,
                                                           filenames=["output/%s%d.json" % (input_filename_base, x)
                                                                      for x in range(start_year, end_year + 1)]
                                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.output_filename = output_filename

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year + 1), self.filenames):
            # from the 30 day reporting interval, imagine all months have 30 days, except December, which has 35
            days_in_month = [30] * 11 + [35]

            # population size
            pop = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][
                  :12]  # remove final five days: assume final five days have same average as rest of month
            pop_Under15 = [sum(x[:3]) for x in pop]
            pop_15to30 = [x[3] for x in pop]
            pop_30to50 = [x[4] for x in pop]
            pop_50plus = [x[5] for x in pop]

            # new infections
            d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:13]
            d[11] = [sum(x) for x in zip(d[11], d[12])]  # add final five days to last month
            del d[-1]  # remove final five days
            new_infections_Under15 = [sum(x[:3]) for x in d]
            new_infections_15to30 = [x[3] for x in d]
            new_infections_30to50 = [x[4] for x in d]
            new_infections_50plus = [x[5] for x in d]

            # PfPR
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # use weighted average for combined age groups
            pfpr_Under15 = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2])) / (
                    pop[yy][0] + pop[yy][1] + pop[yy][2]) for yy in range(12)]
            pfpr_15to30 = [x[3] for x in d]
            pfpr_30to50 = [x[4] for x in d]
            pfpr_50plus = [x[5] for x in d]

            # clinical cases
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # adjust the per-person annualized number (the reported value) to get the total number of clinical cases in that age group in a month
            clinical_cases_Under15 = [
                ((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2])) * days_in_month[yy] / 365
                for yy in range(12)]
            clinical_cases_15to30 = [d[yy][3] * pop[yy][3] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_30to50 = [d[yy][4] * pop[yy][4] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_50plus = [d[yy][5] * pop[yy][5] * days_in_month[yy] / 365 for yy in range(12)]

            # severe cases
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # adjust the per-person annualized number (the reported value) to get the total number of severe cases in that age group in a month
            severe_cases_Under15 = [
                ((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2])) * days_in_month[yy] / 365
                for yy in range(12)]
            severe_cases_15to30 = [d[yy][3] * pop[yy][3] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_30to50 = [d[yy][4] * pop[yy][4] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_50plus = [d[yy][5] * pop[yy][5] * days_in_month[yy] / 365 for yy in range(12)]

            # order is [under 15, 15-30, 30-50, over 50]
            simdata = pd.DataFrame({'month': list(range(1, 13)) * 4,  # cycle through months for each age range
                                    'AgeGroup': np.repeat(['Under15', '15to30', '30to50', '50plus'], 12),
                                    'Pop': (pop_Under15 + pop_15to30 + pop_30to50 + pop_50plus),
                                    'New Infections': (
                                            new_infections_Under15 + new_infections_15to30 + new_infections_30to50 + new_infections_50plus),
                                    'PfPR': (pfpr_Under15 + pfpr_15to30 + pfpr_30to50 + pfpr_50plus),
                                    'Clinical cases': (
                                            clinical_cases_Under15 + clinical_cases_15to30 + clinical_cases_30to50 + clinical_cases_50plus),
                                    'Severe cases': (
                                            severe_cases_Under15 + severe_cases_15to30 + severe_cases_30to50 + severe_cases_50plus)
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, self.output_filename)), index=False)


class MonthlyNewInfectionsAnalyzer_withU5(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026,
                 input_filename_base='MalariaSummaryReport_Monthly',
                 output_filename='newInfections_PfPR_cases_monthly_byAgeGroup_withU5.csv'):

        super(MonthlyNewInfectionsAnalyzer_withU5, self).__init__(working_dir=working_dir,
                                                                  filenames=[
                                                                      "output/%s%d.json" % (input_filename_base, x)
                                                                      for x in range(start_year, end_year + 1)]
                                                                  )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.output_filename = output_filename

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year + 1), self.filenames):
            # from the 30 day reporting interval, imagine all months have 30 days, except December, which has 35
            days_in_month = [30] * 11 + [35]

            # population size
            pop = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][
                  :12]  # remove final five days: assume final five days have same average as rest of month
            pop_Under5 = [sum(x[:2]) for x in pop]
            pop_5to15 = [x[2] for x in pop]
            pop_15to30 = [x[3] for x in pop]
            pop_30to50 = [x[4] for x in pop]
            pop_50plus = [x[5] for x in pop]
            pop_allAges = [sum(x[:6]) for x in pop]

            # new infections
            d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:13]
            d[11] = [sum(x) for x in zip(d[11], d[12])]  # add final five days to last month
            del d[-1]  # remove final five days
            new_infections_Under5 = [sum(x[:2]) for x in d]
            new_infections_5to15 = [x[2] for x in d]
            new_infections_15to30 = [x[3] for x in d]
            new_infections_30to50 = [x[4] for x in d]
            new_infections_50plus = [x[5] for x in d]
            new_infections_allAges = [sum(x[:6]) for x in d]

            # PfPR
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # use weighted average for combined age groups
            pfpr_Under5 = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1])) / (pop[yy][0] + pop[yy][1]) for yy in
                           range(12)]
            pfpr_5to15 = [x[2] for x in d]
            pfpr_15to30 = [x[3] for x in d]
            pfpr_30to50 = [x[4] for x in d]
            pfpr_50plus = [x[5] for x in d]
            pfpr_allAges = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2]) +
                             (d[yy][3] * pop[yy][3]) + (d[yy][4] * pop[yy][4]) + (d[yy][5] * pop[yy][5])) /
                            (pop[yy][0] + pop[yy][1] + pop[yy][2] + pop[yy][3] + pop[yy][4] + pop[yy][5]) for yy in
                            range(12)]

            # clinical cases
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # adjust the per-person annualized number (the reported value) to get the total number of clinical cases in that age group in a month
            clinical_cases_Under5 = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1])) * days_in_month[yy] / 365 for
                                     yy in range(12)]
            clinical_cases_5to15 = [d[yy][2] * pop[yy][2] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_15to30 = [d[yy][3] * pop[yy][3] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_30to50 = [d[yy][4] * pop[yy][4] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_50plus = [d[yy][5] * pop[yy][5] * days_in_month[yy] / 365 for yy in range(12)]
            clinical_cases_allAges = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2]) +
                                       (d[yy][3] * pop[yy][3]) + (d[yy][4] * pop[yy][4]) + (d[yy][5] * pop[yy][5])) *
                                      days_in_month[yy] / 365 for yy in range(12)]

            # severe cases
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][
                :12]  # remove final five days: assume final five days have same average as rest of month
            # adjust the per-person annualized number (the reported value) to get the total number of severe cases in that age group in a month
            severe_cases_Under5 = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1])) * days_in_month[yy] / 365 for yy
                                   in range(12)]
            severe_cases_5to15 = [d[yy][2] * pop[yy][2] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_15to30 = [d[yy][3] * pop[yy][3] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_30to50 = [d[yy][4] * pop[yy][4] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_50plus = [d[yy][5] * pop[yy][5] * days_in_month[yy] / 365 for yy in range(12)]
            severe_cases_allAges = [((d[yy][0] * pop[yy][0]) + (d[yy][1] * pop[yy][1]) + (d[yy][2] * pop[yy][2]) +
                                     (d[yy][3] * pop[yy][3]) + (d[yy][4] * pop[yy][4]) + (d[yy][5] * pop[yy][5])) *
                                    days_in_month[yy] / 365 for yy in range(12)]

            # order is [under 15, 15-30, 30-50, over 50]
            simdata = pd.DataFrame({'month': list(range(1, 13)) * 6,  # cycle through months for each age range
                                    'AgeGroup': np.repeat(['Under5', '5to15', '15to30', '30to50', '50plus', 'allAges'],
                                                          12),
                                    'Pop': (
                                            pop_Under5 + pop_5to15 + pop_15to30 + pop_30to50 + pop_50plus + pop_allAges),
                                    'New Infections': (
                                            new_infections_Under5 + new_infections_5to15 + new_infections_15to30 + new_infections_30to50 + new_infections_50plus + new_infections_allAges),
                                    'PfPR': (
                                            pfpr_Under5 + pfpr_5to15 + pfpr_15to30 + pfpr_30to50 + pfpr_50plus + pfpr_allAges),
                                    'Clinical cases': (
                                            clinical_cases_Under5 + clinical_cases_5to15 + clinical_cases_15to30 + clinical_cases_30to50 + clinical_cases_50plus + clinical_cases_allAges),
                                    'Severe cases': (
                                            severe_cases_Under5 + severe_cases_5to15 + severe_cases_15to30 + severe_cases_30to50 + severe_cases_50plus + severe_cases_allAges)
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

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, self.output_filename)), index=False)


class MonthlyNewInfectionsAnalyzerByAge(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026,
                 input_filename_base='MalariaSummaryReport_Monthly',
                 output_filename='newInfections_PfPR_cases_monthly_byAgeGroup.csv'):

        super(MonthlyNewInfectionsAnalyzerByAge, self).__init__(working_dir=working_dir,
                                                                filenames=["output/%s%d.json" % (input_filename_base, x)
                                                                           for x in range(start_year, end_year + 1)]
                                                                )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.output_filename = output_filename

    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def map(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year + 1), self.filenames):

            # from the 30 day reporting interval, imagine all months have 30 days, except December, which has 35
            days_in_month = [30] * 11 + [35]

            # iterate through age bins, extracting the monthly values of each metric and then appending into data frame
            simdata_allAges = pd.DataFrame()

            for aa in range(len(data[fname]['Metadata']['Age Bins'])):
                # population size
                pop = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][
                      :12]  # remove final five days: assume final five days have same average as rest of month
                pop_monthly = [x[aa] for x in pop]

                # new infections
                d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:13]
                d[11] = [sum(x) for x in zip(d[11], d[12])]  # add final five days to last month
                del d[-1]  # remove final five days
                new_infections_monthly = [x[aa] for x in d]

                # PfPR
                d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][
                    :12]  # remove final five days: assume final five days have same average as rest of month
                pfpr_monthly = [x[aa] for x in d]

                # clinical cases
                d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][
                    :12]  # remove final five days: assume final five days have same average as rest of month
                # adjust the per-person annualized number (the reported value) to get the total number of clinical cases in that age group in a month
                clinical_cases_monthly = [d[yy][aa] * pop[yy][aa] * days_in_month[yy] / 365 for yy in range(12)]

                # severe cases
                d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][
                    :12]  # remove final five days: assume final five days have same average as rest of month
                # adjust the per-person annualized number (the reported value) to get the total number of severe cases in that age group in a month
                severe_cases_monthly = [d[yy][aa] * pop[yy][aa] * days_in_month[yy] / 365 for yy in range(12)]

                # order is [under 15, 15-30, 30-50, over 50]
                simdata = pd.DataFrame({'month': list(range(1, 13)),
                                        'AgeGroup': np.repeat([data[fname]['Metadata']['Age Bins'][aa]], 12),
                                        'Pop': pop_monthly,
                                        'New Infections': new_infections_monthly,
                                        'PfPR': pfpr_monthly,
                                        'Clinical cases': clinical_cases_monthly,
                                        'Severe cases': severe_cases_monthly
                                        })
                simdata['year'] = year
                simdata_allAges = pd.concat([simdata_allAges, simdata])
            adf = pd.concat([adf, simdata_allAges])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf

    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        output_dir = os.path.join(self.working_dir, self.expt_name)
        os.makedirs(output_dir, exist_ok=True)

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, self.output_filename)), index=False)


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

    data_path, project_path = load_box_paths(country_name='Nigeria')

    working_dir = os.path.join(project_path, 'simulation_output', 'simulations_to_present')
    start_year = 2010  # simulation starts in January of this year
    end_year = 2024  # simulation ends in December of this year
    # start_year = 2021  # simulation starts in January of this year
    # end_year = 2030  # simulation ends in December of this year

    expt_ids = {
        'NGA25_toPresent_allInter': 'af6003bc-9778-f011-9f17-b88303912b51'
        # 'NGA_toPresent_allInter': '0a297502-088c-ed11-aa00-b88303911bc1',
    }

    for expname, expid in expt_ids.items():
        print('running expt %s' % expname)
        # report_count_channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment',
        #                           'Bednet_Got_New_One', 'Bednet_Using',  # 'Received_Self_Medication',
        #                          'Received_Campaign_Drugs', 'Received_IRS'
        #                          ]
        report_count_channels = None

        if 'no_IRS_SMC_ITN_CM' in expname:
            cur_monthlyTreatedCasesAnalyzer = monthlyTreatedCasesAnalyzer(expt_name=expname,
                                                                          channels=['Received_NMF_Treatment'],
                                                                          sweep_variables=["Run_Number",
                                                                                           "admin_name"],
                                                                          working_dir=working_dir,
                                                                          start_year=start_year,
                                                                          end_year=end_year)
        else:
            cur_monthlyTreatedCasesAnalyzer = monthlyTreatedCasesAnalyzer(expt_name=expname,
                                                                          channels=report_count_channels,
                                                                          sweep_variables=["Run_Number",
                                                                                           "admin_name"],
                                                                          working_dir=working_dir,
                                                                          start_year=start_year,
                                                                          end_year=end_year)

        analyzers = [
            monthlyU5PfPRAnalyzer(expt_name=expname,
                                  sweep_variables=["Run_Number", "admin_name"],
                                  working_dir=working_dir,
                                  start_year=start_year,
                                  end_year=end_year),
            # # ==== <- remove U1 for 2010-2020 if no IPTi
            # monthlyU1PfPRAnalyzer(expt_name=expname,
            #                       sweep_variables=["Run_Number", "admin_name"],
            #                       working_dir=working_dir,
            #                       start_year=start_year,
            #                       end_year=end_year),
            # # =====
            cur_monthlyTreatedCasesAnalyzer,
            monthlyEventAnalyzer(expt_name=expname,
                                 channels=report_count_channels,
                                 sweep_variables=["Run_Number", "admin_name"],
                                 working_dir=working_dir,
                                 start_year=start_year,
                                 end_year=end_year),
            monthlySevereTreatedByAgeAnalyzer(expt_name=expname,
                                              sweep_variables=["Run_Number", "admin_name"],
                                              working_dir=working_dir,
                                              start_year=start_year,
                                              end_year=end_year,
                                              agebins=[5, 120]),
            MonthlyNewInfectionsAnalyzer(expt_name=expname,
                                         sweep_variables=["Run_Number", "admin_name"],
                                         working_dir=working_dir,
                                         start_year=start_year,
                                         end_year=end_year,
                                         input_filename_base='MalariaSummaryReport_Monthly',
                                         output_filename='newInfections_PfPR_cases_monthly_byAgeGroup.csv'),
            MonthlyNewInfectionsAnalyzer_withU5(expt_name=expname,
                                                sweep_variables=["Run_Number", "admin_name"],
                                                working_dir=working_dir,
                                                start_year=start_year,
                                                end_year=end_year,
                                                input_filename_base='MalariaSummaryReport_Monthly',
                                                output_filename='newInfections_PfPR_cases_monthly_byAgeGroup_withU5.csv'),
            MonthlyNewInfectionsAnalyzerByAge(expt_name=expname,
                                              sweep_variables=["Run_Number", "admin_name"],
                                              working_dir=working_dir,
                                              start_year=start_year,
                                              end_year=end_year,
                                              input_filename_base='MalariaSummaryReport_Monthly',
                                              output_filename='newInfections_PfPR_cases_monthly_byAgeGroup.csv')

        ]
        am = AnalyzeManager(platform=platform, ids=[(expid, ItemType.EXPERIMENT)], analyzers=analyzers,
                            analyze_failed_items=True)
        am.analyze()

