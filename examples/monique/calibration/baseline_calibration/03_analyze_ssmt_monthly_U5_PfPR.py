# running the Platform analyzer will create output files on COMPS in the Work Items section.
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from snt.analyzers.analyze_monthly_pfpr_u5 import MonthlyPfPRU5Analyzer
import tempfile
from idmtools.core import ItemType

experiments = {'PfPR_sweep_main_example': 'e369e323-5153-ee11-aa0a-b88303911bc1'}

start_year = 2010
end_year = 2017

working_dir = "."
print(f"\n{working_dir}\n")


if __name__ == "__main__":
    platform = Platform('Calculon')
    # platform = Platform('IDMCLOUD')

    analyzers = [MonthlyPfPRU5Analyzer]

    sweep_variables = ["Run_Number",
                       "Habitat_Multiplier",
                       "admin_name",
                       ]

    for expt_name, exp_id in experiments.items():
        # if not os.path.exists(os.path.join(working_dir, expt_name)):
        #     os.mkdir(os.path.join(working_dir, expt_name))
        # working_dir = os.path.join(working_dir, expt_name)

        wi_name = "ssmt_analyzer_%s" % expt_name

        args_each = {'expt_name': expt_name,
                     'sweep_variables': sweep_variables,
                     'working_dir': working_dir,
                     'start_year': start_year,
                     'end_year': end_year
                     }
        analysis = PlatformAnalysis(platform=platform,
                                    experiment_ids=[exp_id],
                                    analyzers=analyzers,
                                    analyzers_args=[args_each] * len(analyzers),
                                    analysis_name=wi_name)

        analysis.analyze()
        wi = analysis.get_work_item()
        dirpath = tempfile.mkdtemp()
        out_filenames = [f"monthly_U5_PfPR.csv"]
        platform.get_files_by_id(wi.uid, ItemType.WORKFLOW_ITEM, out_filenames, dirpath)
