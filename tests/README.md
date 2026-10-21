# Tests

```
powerflow_tests          # all unit and integration tests, one binary
ctest --test-dir build   # the unit tests plus a set of CLI smoke tests
```

The suite needs no framework: `test_framework.hpp` is a ~120-line registry with
`TEST`, `CHECK`, `CHECK_EQ`, `CHECK_NEAR` and `CHECK_THROWS_AS`. `main` runs
every registered case and returns a non-zero exit code when anything fails, so
CTest and CI pick failures up directly.

## Layout

| file | covers |
|------|--------|
| `test_dense.cpp` | LU factorisation with partial pivoting, singularity detection, `matvec` |
| `test_case_file.cpp` | the MATPOWER reader, including every error path |
| `test_network.cpp` | Y-bus assembly, the transformer branch model, bus classification, input validation |
| `test_power_flow.cpp` | the solvers, the analytic Jacobian, loss conservation |
| `test_report.cpp` | tables, CSV, JSON and the SVG renderer |

## Where the expected numbers come from

The tests would only prove that the code agrees with itself, so the loss,
voltage and branch-flow expectations are **externally anchored**:

1. **Losses and bus voltages.** A separate implementation of the same equations
   was written in Python with numpy (`Newton-Raphson`, finite-difference
   Jacobian, different internal structure) and run on the same case files. Its
   answers are what `kReferences` in `test_power_flow.cpp` holds. Independently
   of that, the IEEE 14-bus result was checked against the published solution:
   the solved voltages reproduce the values stored in the IEEE Common Data
   Format to 1.3e-3 pu and 0.017°, and the active loss is 13.3933 MW against
   13.39 MW in the literature.

2. **The Jacobian.** `JacobianMatchesFiniteDifferences` compares every entry of
   the analytic Jacobian against a central difference of the mismatch function
   at a realistic operating point. This is the test that matters most, because
   the analytic Jacobian is the easiest thing in the project to get subtly
   wrong, and a wrong Jacobian still often *converges* — just not quadratically.

3. **Conservation.** `BusInjectionEqualsBranchFlowsPlusShunts` checks that the
   sum of the bus injections equals the sum of the branch losses plus the power
   absorbed by the bus shunts. That single identity ties the Y-bus assembly, the
   branch model, the flows and the injections together; it holds to 1e-12.

4. **A second solver.** `GaussSeidelAgreesWithNewton` solves the same case with
   a completely different iteration and requires the same voltages to 1e-6 pu.
   It also asserts that Newton needs strictly fewer iterations, which documents
   the point of having both.

Reference losses (MW / MVAr) at a 1e-12 pu tolerance:

| case | P loss | Q loss |
|------|-------:|-------:|
| IEEE 9   |    4.6410 |  -92.1601 |
| IEEE 14  |   13.3933 |   30.1224 |
| IEEE 30  |    2.4438 |   -6.5627 |
| IEEE 57  |   27.8638 |    6.3280 |
| IEEE 118 |  132.8629 | -557.9474 |

Negative reactive loss is expected on the 9- and 118-bus systems: their lines
carry enough shunt capacitance that the network generates more reactive power
than it consumes.
