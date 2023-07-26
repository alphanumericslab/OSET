import unittest

from lp_filter_zero_phase_unit_test import lp_filter_zero_phase_unit_test
from peak_detection_amp_threshold_unit_test import peak_detection_amp_threshold_unit_test
from peak_detection_local_search_unit_test import peak_detection_local_search_unit_test
from peak_detection_modified_pan_tompkins_unit_test import peak_detection_modified_pan_tompkins_unit_test
from peak_detection_pan_tompkins_unit_test import peak_detection_pan_tompkins_unit_test
from peak_detection_simple_unit_test import peak_detection_simple_unit_test


class TestMyFunctions(unittest.TestCase):
    def test_lp_filter_zero_phase(self):
        self.assertTrue(lp_filter_zero_phase_unit_test())

    def test_peak_detection_amp_threshold(self):
        self.assertTrue(peak_detection_amp_threshold_unit_test())

    def test_peak_detection_local_search(self):
        self.assertTrue(peak_detection_local_search_unit_test())

    def test_peak_detection_modified_pan_tompkins(self):
        self.assertTrue(peak_detection_modified_pan_tompkins_unit_test())

    def test_peak_detection_pan_tompkins(self):
        self.assertTrue(peak_detection_pan_tompkins_unit_test())

    def test_peak_detection_simple(self):
        self.assertTrue(peak_detection_simple_unit_test())


if __name__ == '__main__':
    unittest.main()
