# Test-system data

The five `.m` files in this directory are **MATPOWER-style case files** holding
the public IEEE benchmark networks. They are inputs to `powerflow`, not code.

| file | system | buses | generators | branches |
|------|--------|------:|-----------:|---------:|
| `case9.m`   | WSCC 9-bus      |   9 |  3 |   9 |
| `case14.m`  | IEEE 14-bus     |  14 |  5 |  20 |
| `case30.m`  | IEEE 30-bus     |  30 |  6 |  41 |
| `case57.m`  | IEEE 57-bus     |  57 |  7 |  80 |
| `case118.m` | IEEE 118-bus    | 118 | 54 | 186 |

## Format

The subset of the MATPOWER case format that a power flow needs:

```
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus    = [ ... ];   % 13 columns, at least
mpc.gen    = [ ... ];   % 10 columns, at least
mpc.branch = [ ... ];   % 13 columns, at least
```

Other sections that real MATPOWER files carry (`mpc.gencost`, `mpc.bus_name`,
`mpc.gentype`, ...) are ignored by the parser, so a stock MATPOWER case file
loads unchanged. The files here were slimmed down to keep the repository small;
the numbers are untouched.

Column meanings are documented in `include/powerflow/types.hpp` and in
`src/case_file.cpp`.

## Provenance

The numerical values are the IEEE test systems as published by the
**University of Washington Power Systems Test Case Archive**
(<https://labs.ece.uw.edu/pstca/>) in IEEE Common Data Format, converted to the
MATPOWER layout by the MATPOWER project
(<https://github.com/MATPOWER/matpower>, BSD-3-Clause). The conversion script
lives outside the repository; only the resulting numbers were carried over, so
no third-party source text is redistributed here.

## Can I add my own case?

Yes. The easiest route is to export from MATPOWER/PYPOWER or to write the three
matrices by hand. Requirements the solver enforces:

* exactly one bus of type 3 (slack) among the connected buses,
* every branch endpoint refers to an existing bus id,
* no branch with `r = x = 0` in service,
* the network does not have to be fully meshed, but every bus that should take
  part in the solution needs at least one in-service branch (otherwise it is
  reported as isolated).
