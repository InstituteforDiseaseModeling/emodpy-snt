import logging
from abc import ABCMeta
from idmtools_calibra.calib_site import CalibSite
from snt.calibration.helpers_seasonality_calibration import get_cases
from snt.calibration.ChannelByMultiYearSeasonCohortInsetAnalyzer import ChannelByMultiYearSeasonCohortInsetAnalyzer

logger = logging.getLogger(__name__)


class SeasonalityCalibSite(CalibSite):
    __metaclass__ = ABCMeta

    def __init__(self, hfca=None, project_path=None, **kwargs):
        self.metadata = {
            'HFCA': hfca
        }
        if 'throwaway' in kwargs:
            self.throwaway = kwargs['throwaway']
        else:
            self.throwaway = 0
        self.name = hfca
        self.project_path = project_path

        super(SeasonalityCalibSite, self).__init__(hfca)

    def get_reference_data(self, reference_type):
        super(SeasonalityCalibSite, self).get_reference_data(reference_type)

        # Load the case CSV
        reference_data = get_cases(self.metadata['HFCA'], self.project_path)

        return reference_data

    def get_analyzers(self):

        return [
            ChannelByMultiYearSeasonCohortInsetAnalyzer(site=self, throwaway=self.throwaway)]
