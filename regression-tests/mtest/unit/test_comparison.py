#!/usr/bin/env python
"""
Unit test suite for the betr regression test manager

"""
from __future__ import print_function

import logging
import os
import os.path
import sys
import unittest

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


from rtest_betr import Comparison


class Comparison_suite(unittest.TestCase):
    """
    """
    _LOG_FILENAME = 'dummy.testlog'

    def setUp(self):
        """
        """
        logging.basicConfig(filename=self._LOG_FILENAME,
                            filemode='w',
                            level=logging.INFO,
                            format='%(message)s')
        logging.info('mtest {0} unit test log.'.format(__name__))

    def tearDown(self):
        """
        """
        logging.shutdown()
        if os.path.isfile(self._LOG_FILENAME):
            os.remove(self._LOG_FILENAME)

    # ------------------------------------------------------

    def test_comparison_get_section_category_valid(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)

        section = {'category': comparison._tolerances.CONC,
                   'foo': 'bar baz'}
        expected = comparison._tolerances.CONC
        received = comparison._get_section_category('some test', section)
        self.assertEqual(expected, received)

    def test_comparison_get_section_category_missing_category(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)

        section = {'pet': 'dog cat',
                   'foo': 'bar baz'}

        self.assertRaises(RuntimeError,
                          comparison._get_section_category,
                          'some test', section)

    def test_comparison_get_section_category_invalid_category(self):
        """Test that we can extract the category correctly from a test
        description section.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        
        section = {'category': 'dog cat',
                   'foo': 'bar baz'}

        self.assertRaises(RuntimeError,
                          comparison._get_section_category,
                          'some test', section)

    # ------------------------------------------------------
    def test_comparison_float_absolute_pass(self):
        """Test that comparison of float with absolute tolerance passes when
        it is less that tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-16 absolute')
        section = 'Ca'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.1e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertTrue(received)

    def test_comparison_float_absolute_fail(self):
        """Test that comparison of float with absolute tolerance fails when it
        is greater than tolerance.

        """
        conf = {}
        comparison = Comparison('unittest', conf)
        category = comparison._tolerances.CONC
        comparison.update_from_name(category, '1.0e-18 absolute')
        section = 'Ca'
        key = 'cell 1'
        a_data = '1.0e-16'
        b_data = '1.1e-16'

        received = comparison._compare_float_values_with_tolerance(
            category, section, key, a_data, b_data)
        self.assertFalse(received)

if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()
