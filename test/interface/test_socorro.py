import unittest

import numpy as np
from phonopy.interface.phonopy_yaml import get_unitcell_from_phonopy_yaml
from phonopy.interface.socorro import read_socorro
import os

data_dir = os.path.dirname(os.path.abspath(__file__))

#
# need NaCl-socorro-crystal, which should be read and then output as NaCl-abinit-pwscf.yaml
#

class TestSocorro(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_socorro(self):
        cell = read_socorro(os.path.join(data_dir, "NaCl-socorro.in"))
        filename = os.path.join(data_dir, "NaCl-abinit-pwscf.yaml")
        cell_ref = get_unitcell_from_phonopy_yaml(filename)
        print cell_ref.get_scaled_positions()
        self.assertTrue(
            (np.abs(cell.get_cell() - cell_ref.get_cell()) < 1e-5).all())
        diff_pos = (cell.get_scaled_positions()
                    - cell_ref.get_scaled_positions())
        diff_pos -= np.rint(diff_pos)
        self.assertTrue((np.abs(diff_pos) < 1e-5).all())
        for s, s_r in zip(cell.get_chemical_symbols(),
                          cell_ref.get_chemical_symbols()):
            self.assertTrue(s == s_r)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSocorro)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
