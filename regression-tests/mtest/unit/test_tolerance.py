#!/usr/bin/env python
"""
Unit test suite for the betr regression test manager

"""
from __future__ import print_function

import logging
import sys
import unittest

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser


from rtest_betr import Tolerances


class Tolerances_suite(unittest.TestCase):
    """
    """

    def setUp(self):
        """
        """
        log_filename = 'dummy.testlog'
        logging.basicConfig(filename=log_filename,
                            filemode='w',
                            level=logging.INFO,
                            format='%(message)s')

    def tearDown(self):
        """
        """
        logging.shutdown()

    # ------------------------------------------------------

    def test_tolerances_default(self):
        """Test basic initialization of a tolerance object and return of a
        default value.

        """
        tolerances = Tolerances()
        expected = Tolerances._DEFAULT_EPSILON
        received = tolerances.get(Tolerances.CONC, 'value')

        self.assertEqual(expected, received)


if __name__ == '__main__':
    # unittest.main(buffer=True)
    unittest.main()
