13.May.2025

At some later date, I rewrote all of this functionality from scratch so that I could use bit-fiddling on occupation strings
on the C backend for efficiency.  I was pretty happy with that, and wanted to deprecate this earlier version.
First I looked to see where it is still in use in live code (not in an attic or old repo), and it turns out it was
already out of use.  Below is every reference to fermion_field in a .py file in a subdirectory of Qode before the move
was made.

Development/old-repos/Qode-BeStates/qode/__init__.py:from . import fermion_field
Development/old-repos/Qode-optimizer/qode/__init__.py:from . import fermion_field
Development/old-repos/Qode/qode/__init__.py:from . import fermion_field
qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings,0)
qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings)
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.state import state
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.nested_operator2 import nested_operator as operator
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.reorder2 import *
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field          import occ_strings
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.state    import configuration, state, dot
qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.reorder2 import *
qode/fermion_field/attic/hopefully_deprecated/update_occ_strings/test.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode-BeStates/Applications/GPU-pilot/attic/debug_Hbuild.py:from qode.fermion_field.state import *
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.state import state
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field import nested_operator as operator
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test1.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field.state import *
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.state   import configuration, state, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.nested_operator import multiply as op_mult
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.state   import configuration, state, dot, op_string, create, annihilate
Development/old-repos/Qode-BeStates/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-optimizer/Applications/GPU-pilot/attic/debug_Hbuild.py:from qode.fermion_field.state import *
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.state import state
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field import nested_operator as operator
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test1.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field.state import *
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.state   import configuration, state, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.nested_operator import multiply as op_mult
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.state   import configuration, state, dot, op_string, create, annihilate
Development/old-repos/Qode-optimizer/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode/Applications/GPU-pilot/attic/debug_Hbuild.py:from qode.fermion_field.state import *
Development/old-repos/Qode/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/He_FCI.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.state import state
Development/old-repos/Qode/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field import nested_operator as operator
Development/old-repos/Qode/Applications/component_tests/test_fermions/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode/Applications/component_tests/test_fermions/test1.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/test2.py:from qode.fermion_field.state import *
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.state   import configuration, state, dot
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_mult.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_op_multiply.py:from qode.fermion_field.nested_operator import multiply as op_mult
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field         import occ_strings
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.state   import configuration, state, dot, op_string, create, annihilate
Development/old-repos/Qode/Applications/component_tests/test_fermions/test_reorder.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field       import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/elec_cisd_main.py:from qode.fermion_field import occ_strings, state
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate, dot
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/hamiltonian_simple.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field.state import state, dot, resolvent
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:	configurations = qode.fermion_field.occ_strings.CISD(occ_orbs,vrt_orbs)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/t_operator.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import excitation    as Ex
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import deexcitation  as Dx
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import occ_rearrange as Fo
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import vrt_rearrange as Fv
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/local/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/recoupled-osc-ccsd/attic/original/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-BeStates/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings,0)
Development/old-repos/Qode-BeStates/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field       import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/elec_cisd_main.py:from qode.fermion_field import occ_strings, state
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate, dot
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/hamiltonian_simple.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field.state import state, dot, resolvent
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:	configurations = qode.fermion_field.occ_strings.CISD(occ_orbs,vrt_orbs)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/t_operator.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import excitation    as Ex
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import deexcitation  as Dx
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import occ_rearrange as Fo
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import vrt_rearrange as Fv
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/local/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/recoupled-osc-ccsd/attic/original/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-optimizer/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings,0)
Development/old-repos/Qode-optimizer/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings)
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/build_H_from_class_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field       import occ_strings
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/build_H_from_function_and_diag.py:from qode.fermion_field.state import state, dot
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/elec_cisd_main.py:from qode.fermion_field import occ_strings, state
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate, dot
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/hamiltonian.py:from qode.fermion_field.reorder import *
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/hamiltonian_simple.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field import occ_strings
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/non_linear_opt_w_state_arithmetic.py:from qode.fermion_field.state import state, dot, resolvent
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:from qode.fermion_field.state import state, op_string, create, annihilate
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/simpleHamiltonian.py:	configurations = qode.fermion_field.occ_strings.CISD(occ_orbs,vrt_orbs)
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/t_operator.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import excitation    as Ex
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import deexcitation  as Dx
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import occ_rearrange as Fo
Development/old-repos/Qode/Applications/component_tests/ccsd/fast/hamiltonian.py:from qode.fermion_field.reorder import vrt_rearrange as Fv
Development/old-repos/Qode/Applications/component_tests/ccsd/local/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/ccsd.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/general-XRCC/attic/FCI_tests/hard_way.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core2.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core3.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/general-XRCC/attic/core_test/rm_core4.py:from qode.fermion_field.state import create, annihilate, op_string
Development/old-repos/Qode/Applications/recoupled-osc-ccsd/attic/original/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings,0)
Development/old-repos/Qode/qode/many_body/coupled_cluster/attic/operator_commutation.py:	return_state = qode.fermion_field.state.state(occ_strings)
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/mol_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/mol_fci/excite_strings.py:from qode.fermion_field.occ_strings import all_occ_strings as excited_molecules		# a coincident usage of same logic
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-BeStates/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.state       import configuration, state, create, annihilate, op_string
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/bch.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/build_H_into_blocks.py:# from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/hamiltonian.py:from qode.fermion_field import state, config
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/main.py:from qode.fermion_field import state
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode                                       import fermion_field
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/t_operator.py:from qode.fermion_field import state, config
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:from qode                           import fermion_field
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:	fluctuations = fermion_field.OV_partitioning(orbitals.occ, orbitals.vrt)
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/util.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/build_Omega.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/check_decomp.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/decompose_commutator.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/diff_ops.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/sum.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/nested/debug-util-elec/view_op.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-BeStates/Applications/recoupled-osc-ccsd/attic/validate/old/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.state import state
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.nested_operator2 import nested_operator as operator
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field          import occ_strings
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.state    import configuration, state, dot
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode-BeStates/qode/fermion_field/attic/hopefully_deprecated/update_occ_strings/test.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/mol_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/mol_fci/excite_strings.py:from qode.fermion_field.occ_strings import all_occ_strings as excited_molecules		# a coincident usage of same logic
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode-optimizer/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.state       import configuration, state, create, annihilate, op_string
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/bch.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/build_H_into_blocks.py:# from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/hamiltonian.py:from qode.fermion_field import state, config
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/main.py:from qode.fermion_field import state
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode                                       import fermion_field
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/t_operator.py:from qode.fermion_field import state, config
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:from qode                           import fermion_field
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:	fluctuations = fermion_field.OV_partitioning(orbitals.occ, orbitals.vrt)
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/util.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/build_Omega.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/check_decomp.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/decompose_commutator.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/diff_ops.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/sum.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/nested/debug-util-elec/view_op.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode-optimizer/Applications/recoupled-osc-ccsd/attic/validate/old/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.state import state
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.nested_operator2 import nested_operator as operator
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field          import occ_strings
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.state    import configuration, state, dot
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode-optimizer/qode/fermion_field/attic/hopefully_deprecated/update_occ_strings/test.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_2_noCT_NO.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/elec_fci/H2_counterpoise.py:from qode.fermion_field.state       import state, configuration
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/mol_fci/CI_Hamiltonian.py:from qode.fermion_field.state import *
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/mol_fci/excite_strings.py:from qode.fermion_field.occ_strings import all_occ_strings as excited_molecules		# a coincident usage of same logic
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.occ_strings import all_occ_strings
Development/old-repos/Qode/Applications/Be_n/attic/H2_2_FCI/mol_fci/use_rho.py:from qode.fermion_field.state       import configuration, state, create, annihilate, op_string
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/bch.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/build_H_into_blocks.py:# from qode.fermion_field import state
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/hamiltonian.py:from qode.fermion_field import state, config
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/main.py:from qode.fermion_field import state
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:from qode                                       import fermion_field
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/new_hamiltonian.py:	fluctuations = fermion_field.OV_partitioning(occ_orbs,vrt_orbs)
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/t_operator.py:from qode.fermion_field import state, config
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/hamiltonian.py:from qode.fermion_field.reorder     import reorder_CpAq, reorder_CpCqArAs
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:from qode                           import fermion_field
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/main.py:	fluctuations = fermion_field.OV_partitioning(orbitals.occ, orbitals.vrt)
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/CIS-SCF_based-on-deprecated-main/util.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/build_Omega.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/check_decomp.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/decompose_commutator.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/diff_ops.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/sum.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/component_tests/ccsd/nested/debug-util-elec/view_op.py:from qode.fermion_field import OV_partitioning
Development/old-repos/Qode/Applications/recoupled-osc-ccsd/attic/validate/old/ccsd.py:from qode                               import fermion_field
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.state import state
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.nested_operator2 import nested_operator as operator
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/hamiltonian2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field          import occ_strings
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.state    import configuration, state, dot
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/special-case-exceptions/test_mult2.py:from qode.fermion_field.reorder2 import *
Development/old-repos/Qode/qode/fermion_field/attic/hopefully_deprecated/update_occ_strings/test.py:from qode.fermion_field import config, occ_strings
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ccsd_main.py:from qode.fermion_field import state, config
Development/old-repos/Qode-BeStates/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ho_2mol_cisd.py:from qode.fermion_field import state, config
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ccsd_main.py:from qode.fermion_field import state, config
Development/old-repos/Qode-optimizer/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ho_2mol_cisd.py:from qode.fermion_field import state, config
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ccsd_main.py:from qode.fermion_field import state, config
Development/old-repos/Qode/Applications/component_tests/ccsd/attic/oscillator_ccsd/attic/v1/ho_2mol_cisd.py:from qode.fermion_field import state, config
