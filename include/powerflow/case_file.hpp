// Reader for the MATPOWER .m case format (the parts a power flow needs).
#pragma once

#include <string>

#include "powerflow/types.hpp"

namespace pf {

// Parses `mpc.baseMVA`, `mpc.bus`, `mpc.gen` and `mpc.branch`. Other sections
// such as `mpc.gencost` or `mpc.bus_name` are ignored on purpose.
//
// Throws std::runtime_error with a message that names the offending row when a
// required matrix is missing, a row is too short, or a value is not a number.
CaseData load_case_file(const std::string& path);

// Same as above but reads from an in-memory buffer. `name` is used for display
// and for error messages.
CaseData parse_case_text(const std::string& text, const std::string& name);

}  // namespace pf
