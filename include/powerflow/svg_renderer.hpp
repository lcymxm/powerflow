// Renders a solved case as a self-contained SVG single-line diagram.
//
// The layout is produced by a deterministic Fruchterman-Reingold spring
// embedder, so no bus coordinates are needed and the same image comes out of
// every run and every platform.
#pragma once

#include <string>

#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"

namespace pf {

struct SvgOptions {
    int width = 1200;
    int height = 820;
    std::string title;     // empty: derived from the case
    std::string subtitle;  // empty: derived from the result
};

std::string render_svg(const Network& net, const PowerFlowResult& result,
                       const SvgOptions& options = {});

}  // namespace pf
