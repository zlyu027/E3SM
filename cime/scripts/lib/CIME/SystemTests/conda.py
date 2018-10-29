"""
Implemetation of CIME MCC test: Compares ensemble methods

This does two runs: In the first we run a three member ensemble using the
 MULTI_DRIVER capability, then we run a second single instance case and compare
"""

import sys
import logging

from CIME.utils import expect
from CIME.SystemTests.system_tests_common import SystemTestsCommon

logger = logging.getLogger(__name__)


class CONDA(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """
        super(CONDA, self).__init__(case)

        isconda = 'conda' not in sys.version.lower() or 'continuum' not in sys.version.lower()

        expect(isconda, 'This is not running inside a Anaconda/miniconda environment. '
                        'The system version is: "{}"'.format(sys.version))
