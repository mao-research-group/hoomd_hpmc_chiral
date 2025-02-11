# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Template HOOMD-blue component."""

import hoomd
import hoomd.conftest
from hoomd.hpmc_energy import _hpmc_energy


@hoomd.logging.modify_namespace(('hpmc', 'pair', 'ExamplePair'))
class ExamplePair(hoomd.hpmc.pair.Pair):
    r"""An example pair potential for HPMC.

    Args:
        default_r_cut (float): Default cutoff radius :math:`[\mathrm{length}]`.

    TODO: document the class.

    """

    _cpp_class_name = 'ExamplePairPotential'
    _ext_module = _hpmc_energy

    def __init__(self, default_r_cut=None):
        if default_r_cut is None:
            default_r_cut = float
        else:
            default_r_cut = float(default_r_cut)

        params = hoomd.data.typeparam.TypeParameter(
            'params',
            'particle_types',
            # TODO: Add entries to TypeParameterDict matching those in the C++
            # ParamType.
            hoomd.data.parameterdicts.TypeParameterDict(
                A=float, B=float, r_cut=default_r_cut, len_keys=2
            ),
        )
        self._add_typeparam(params)
