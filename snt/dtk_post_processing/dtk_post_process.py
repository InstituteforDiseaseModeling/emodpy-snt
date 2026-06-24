import re
import os

from dtk_get_burden_functions import monthlyU1PfPRAnalyzer, monthlyU5PfPRAnalyzer, monthlyTreatedCasesAnalyzer, monthlySevereTreatedByAgeAnalyzer, MonthlyNewInfectionsAnalyzer_byAgeGroup_withU1U5

def application( output_path ):
    # determine which years are available for MalariaSummaryReport
    pattern = re.compile(r"^MalariaSummaryReport_Monthly(\d{4})")  # specify the pattern for these files
    all_files = os.listdir(output_path)
    # years numbers from matching filenames
    all_years = [
        int(m.group(1))
        for f in all_files
        if (m := pattern.match(f))
    ]
    start_year, end_year = min(all_years), max(all_years)

    monthlyU1PfPRAnalyzer(output_path, start_year, end_year)
    monthlyU5PfPRAnalyzer(output_path, start_year, end_year)
    monthlyTreatedCasesAnalyzer(output_path, start_year)
    monthlySevereTreatedByAgeAnalyzer(output_path, start_year)
    MonthlyNewInfectionsAnalyzer_byAgeGroup_withU1U5(output_path, start_year, end_year)

if __name__ == "__main__":
    application( "output" )