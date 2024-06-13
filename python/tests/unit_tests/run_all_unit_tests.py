"""This is an automatic script to run all the Python-Matlab Unittests at ones"""

import unittest

from lp_filter_zero_phase_unit_test import lp_filter_zero_phase_unit_test
from peak_det_amp_threshold_unit_test import (
    peak_det_amp_threshold_unit_test,
)
from peak_det_local_search_unit_test import peak_det_local_search_unit_test
from peak_det_matched_filter_robust_unit_test import (
    peak_det_matched_filter_robust_unit_test,
)
from peak_det_matched_filter_unit_test import (
    peak_det_matched_filter_unit_test,
)
from peak_det_modified_pan_tompkins_unit_test import (
    peak_det_modified_pan_tompkins_unit_test,
)
from peak_det_pan_tompkins_unit_test import peak_det_pan_tompkins_unit_test
from peak_det_simple_unit_test import peak_det_simple_unit_test
from baseline_sliding_window_unit_test import baseline_sliding_window_unit_test


class TestMyFunctions(unittest.TestCase):
    def test_lp_filter_zero_phase(self):
        self.assertTrue(lp_filter_zero_phase_unit_test())

    def test_peak_det_amp_threshold(self):
        self.assertTrue(peak_det_amp_threshold_unit_test())

    def test_peak_det_local_search(self):
        self.assertTrue(peak_det_local_search_unit_test())

    def test_peak_det_modified_pan_tompkins(self):
        self.assertTrue(peak_det_modified_pan_tompkins_unit_test())

    def test_peak_det_pan_tompkins(self):
        self.assertTrue(peak_det_pan_tompkins_unit_test())

    def test_peak_det_simple(self):
        self.assertTrue(peak_det_simple_unit_test())

    def test_peak_det_matched_filter_robust(self):
        self.assertTrue(peak_det_matched_filter_robust_unit_test())

    def test_peak_det_matched_filter_unit_test(self):
        self.assertTrue(peak_det_matched_filter_unit_test())

    def test_baseline_sliding_window_unit_test(self):
        self.assertTrue(baseline_sliding_window_unit_test())


if __name__ == "__main__":
    unittest.main()
