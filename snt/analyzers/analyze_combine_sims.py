import os
import re
import datetime
import pandas as pd
import numpy as np
from idmtools.entities import IAnalyzer


class combineSimOutputAnalyzer1(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(combineSimOutputAnalyzer1, self).__init__(working_dir=working_dir,
                                                       filenames=["output/U5_PfPR_ClinicalIncidence_severeTreatment.csv",
                                                                  # "output/U1_PfPR_ClinicalIncidence_severeTreatment.csv",
                                                                  # "output/newInfections_PfPR_cases_monthly_byAgeGroup_withU1U5.csv",
                                                                  # "output/All_Age_monthly_Cases.csv"
                                                                  ]
                                                       )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name


    def map(self, data, simulation):
        adf = data[self.filenames[0]]
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
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, re.sub(r"^output/", "", self.filenames[0]))), index=False)


class combineSimOutputAnalyzer2(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(combineSimOutputAnalyzer2, self).__init__(working_dir=working_dir,
                                                       filenames=[
                                                                  "output/U1_PfPR_ClinicalIncidence_severeTreatment.csv"]
                                                       )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name


    def map(self, data, simulation):
        adf = data[self.filenames[0]]
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
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, re.sub(r"^output/", "", self.filenames[0]))), index=False)


class combineSimOutputAnalyzer3(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(combineSimOutputAnalyzer3, self).__init__(working_dir=working_dir,
                                                        filenames=[
                                                            "output/newInfections_PfPR_cases_monthly_byAgeGroup_withU1U5.csv"]
                                                        )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name

    def map(self, data, simulation):
        adf = data[self.filenames[0]]
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
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, re.sub(r"^output/", "", self.filenames[0]))), index=False)


class combineSimOutputAnalyzer4(IAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(combineSimOutputAnalyzer4, self).__init__(working_dir=working_dir,
                                                        filenames=[
                                                            # "output/U5_PfPR_ClinicalIncidence_severeTreatment.csv",
                                                            # "output/U1_PfPR_ClinicalIncidence_severeTreatment.csv",
                                                            # "output/newInfections_PfPR_cases_monthly_byAgeGroup_withU1U5.csv",
                                                            "output/All_Age_monthly_Cases.csv"]
                                                        )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name

    def map(self, data, simulation):
        adf = data[self.filenames[0]]
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
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, re.sub(r"^output/", "", self.filenames[0]))), index=False)


if __name__ == "__main__":
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform
    from snt.load_paths import load_box_paths

    platform = Platform('Calculon')

    data_path, project_path = load_box_paths(country_name='Nigeria')

    # working_dir = os.path.join(project_path, 'simulation_output', 'simulations_to_present')
    working_dir = os.path.join(project_path, 'simulation_outputs', 'simulations_future')


    expt_ids = {
        'test': '70fbb964-be92-f011-9f19-b88303912b51'
    }

    for expname, expid in expt_ids.items():
        print('running expt %s' % expname)

        analyzers = [
            combineSimOutputAnalyzer1(expt_name=expname,
                                  sweep_variables=["Run_Number", "admin_name"],
                                  working_dir=working_dir),
            combineSimOutputAnalyzer2(expt_name=expname,
                                      sweep_variables=["Run_Number", "admin_name"],
                                      working_dir=working_dir),
            combineSimOutputAnalyzer3(expt_name=expname,
                                      sweep_variables=["Run_Number", "admin_name"],
                                      working_dir=working_dir),
            combineSimOutputAnalyzer4(expt_name=expname,
                                      sweep_variables=["Run_Number", "admin_name"],
                                      working_dir=working_dir),
        ]
        am = AnalyzeManager(platform=platform, ids=[(expid, ItemType.EXPERIMENT)], analyzers=analyzers,
                            analyze_failed_items=True)
        am.analyze()

