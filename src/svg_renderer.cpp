#include "powerflow/svg_renderer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace pf {
namespace {

constexpr const char* kSans =
    "system-ui,-apple-system,'Segoe UI',Roboto,'Helvetica Neue',Arial,sans-serif";
constexpr const char* kMono =
    "ui-monospace,SFMono-Regular,Menlo,Consolas,'Liberation Mono',monospace";

// Dark halo behind label text so it stays readable on top of the branches.
constexpr const char* kHalo = "#0b1220";

struct Rgb {
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
};

struct Point {
    double x = 0.0;
    double y = 0.0;
};

std::string num(double value, int precision = 2) {
    if (std::abs(value) < 5e-13) {
        value = 0.0;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

std::string scientific(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(1) << value;
    return out.str();
}

Rgb mix(Rgb a, Rgb b, double t) {
    const double k = std::clamp(t, 0.0, 1.0);
    return Rgb{a.r + (b.r - a.r) * k, a.g + (b.g - a.g) * k, a.b + (b.b - a.b) * k};
}

std::string hex(Rgb colour) {
    static const char* digits = "0123456789abcdef";
    std::string out = "#";
    for (const double channel : {colour.r, colour.g, colour.b}) {
        const int value = static_cast<int>(std::lround(std::clamp(channel, 0.0, 1.0) * 255.0));
        out.push_back(digits[(value >> 4) & 0xF]);
        out.push_back(digits[value & 0xF]);
    }
    return out;
}

double luminance(Rgb colour) {
    return 0.2126 * colour.r + 0.7152 * colour.g + 0.0722 * colour.b;
}

// Low voltage is red, nominal is green, high voltage is blue.
Rgb voltage_colour(double vm) {
    const Rgb low{0.86, 0.15, 0.15};
    const Rgb nominal{0.13, 0.77, 0.37};
    const Rgb high{0.23, 0.51, 0.96};
    if (vm <= 1.0) {
        return mix(low, nominal, (vm - 0.90) / 0.10);
    }
    return mix(nominal, high, (vm - 1.0) / 0.10);
}

// Cool for a lightly loaded branch, hot for one that is over its rating.
Rgb loading_colour(double ratio, bool rated) {
    if (!rated) {
        return Rgb{0.40, 0.47, 0.58};
    }
    const Rgb cool{0.13, 0.83, 0.93};
    const Rgb comfortable{0.13, 0.77, 0.37};
    const Rgb warm{0.98, 0.80, 0.08};
    const Rgb hot{0.94, 0.27, 0.27};
    if (ratio <= 0.5) {
        return mix(cool, comfortable, ratio / 0.5);
    }
    if (ratio <= 0.85) {
        return mix(comfortable, warm, (ratio - 0.5) / 0.35);
    }
    return mix(warm, hot, (ratio - 0.85) / 0.30);
}

// Deterministic Fruchterman-Reingold layout.
//
// Two adjustments over the textbook version are what make the diagrams usable:
// a weak pull towards the centre, so nodes do not all migrate to the border and
// leave a hole in the middle, and no hard clamping to the box, so the result can
// be scaled uniformly afterwards instead of being squashed per axis.
std::vector<Point> spring_layout(const Network& net, double width, double height,
                                 std::size_t iterations) {
    const std::size_t n = net.bus_count();
    std::vector<Point> position(n);
    if (n == 1) {
        position[0] = Point{width / 2.0, height / 2.0};
        return position;
    }

    std::uint32_t seed = 0x9e3779b9u;
    const auto next_random = [&seed]() {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>(seed >> 8) / static_cast<double>(1u << 24);
    };

    const double centre_x = width / 2.0;
    const double centre_y = height / 2.0;
    const double radius_x = 0.30 * width;
    const double radius_y = 0.30 * height;
    for (std::size_t i = 0; i < n; ++i) {
        const double angle =
            2.0 * 3.14159265358979323846 * static_cast<double>(i) / static_cast<double>(n);
        const double jitter = 0.82 + 0.36 * next_random();
        position[i] = Point{centre_x + radius_x * jitter * std::cos(angle),
                            centre_y + radius_y * jitter * std::sin(angle)};
    }

    struct Edge {
        std::size_t a = 0;
        std::size_t b = 0;
    };
    std::vector<Edge> edges;
    edges.reserve(net.branch_count());
    for (std::size_t k = 0; k < net.branch_count(); ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        edges.push_back(Edge{net.index_of(net.branch(k).from), net.index_of(net.branch(k).to)});
    }

    const double ideal = std::sqrt(width * height / static_cast<double>(n));
    const double gravity = 0.55;  // fraction of the offset applied each iteration
    const double limit_x = 0.75 * width;
    const double limit_y = 0.75 * height;

    std::vector<Point> displacement(n);
    double temperature = 0.14 * std::min(width, height);
    const double cooling =
        iterations > 0 ? std::pow(0.0015 / temperature, 1.0 / static_cast<double>(iterations))
                       : 0.95;

    for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
        for (auto& d : displacement) {
            d = Point{};
        }

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                double dx = position[i].x - position[j].x;
                double dy = position[i].y - position[j].y;
                double distance = std::sqrt(dx * dx + dy * dy);
                if (distance < 1e-3) {
                    dx = 1e-3 * (1.0 + next_random());
                    dy = 1e-3 * (1.0 + next_random());
                    distance = std::sqrt(dx * dx + dy * dy);
                }
                const double force = ideal * ideal / distance;
                const double ux = dx / distance;
                const double uy = dy / distance;
                displacement[i].x += ux * force;
                displacement[i].y += uy * force;
                displacement[j].x -= ux * force;
                displacement[j].y -= uy * force;
            }
        }

        for (const Edge& edge : edges) {
            double dx = position[edge.a].x - position[edge.b].x;
            double dy = position[edge.a].y - position[edge.b].y;
            double distance = std::sqrt(dx * dx + dy * dy);
            if (distance < 1e-3) {
                distance = 1e-3;
            }
            const double force = distance * distance / ideal;
            const double ux = dx / distance;
            const double uy = dy / distance;
            displacement[edge.a].x -= ux * force;
            displacement[edge.a].y -= uy * force;
            displacement[edge.b].x += ux * force;
            displacement[edge.b].y += uy * force;
        }

        const double pull = 0.09 * gravity;
        for (std::size_t i = 0; i < n; ++i) {
            displacement[i].x += (centre_x - position[i].x) * pull;
            displacement[i].y += (centre_y - position[i].y) * pull;

            const double magnitude = std::sqrt(displacement[i].x * displacement[i].x +
                                               displacement[i].y * displacement[i].y);
            if (magnitude > 1e-9) {
                const double step = std::min(magnitude, temperature);
                position[i].x += displacement[i].x / magnitude * step;
                position[i].y += displacement[i].y / magnitude * step;
            }
            // Safety net only: keep the cloud from flying off to infinity.
            position[i].x = std::clamp(position[i].x, centre_x - limit_x, centre_x + limit_x);
            position[i].y = std::clamp(position[i].y, centre_y - limit_y, centre_y + limit_y);
        }

        temperature *= cooling;
    }
    return position;
}

}  // namespace

std::string render_svg(const Network& net, const PowerFlowResult& result,
                       const SvgOptions& options) {
    const int width = std::max(480, options.width);
    const int height = std::max(360, options.height);
    const std::size_t n = net.bus_count();
    const std::size_t nl = net.branch_count();

    const double header_height = 74.0;
    const double footer_height = 132.0;
    const double side_margin = 40.0;
    const double region_x = side_margin;
    const double region_y = header_height;
    const double region_w = static_cast<double>(width) - 2.0 * side_margin;
    const double region_h = static_cast<double>(height) - header_height - footer_height;

    // Run the embedder in a box with the drawing region's aspect ratio, then fit
    // the result uniformly so circles stay round and nothing gets squashed.
    const double layout_w = 1000.0;
    const double layout_h = 1000.0 * region_h / region_w;
    const std::size_t iterations = n > 60 ? 700 : 600;
    const std::vector<Point> layout = spring_layout(net, layout_w, layout_h, iterations);

    double min_x = layout[0].x;
    double max_x = layout[0].x;
    double min_y = layout[0].y;
    double max_y = layout[0].y;
    for (const Point& p : layout) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }
    const double span_x = std::max(max_x - min_x, 1e-6);
    const double span_y = std::max(max_y - min_y, 1e-6);

    const double pad = n > 60 ? 16.0 : 46.0;
    const double avail_w = region_w - 2.0 * pad;
    const double avail_h = region_h - 2.0 * pad;

    // Stretch each axis to fill the frame, but cap the distortion so a very
    // lopsided cloud of buses does not turn into a smear. Positions only, so the
    // bus circles stay perfectly round.
    constexpr double kMaxStretch = 1.6;
    double scale_x = avail_w / span_x;
    double scale_y = avail_h / span_y;
    if (scale_x > scale_y * kMaxStretch) {
        scale_x = scale_y * kMaxStretch;
    }
    if (scale_y > scale_x * kMaxStretch) {
        scale_y = scale_x * kMaxStretch;
    }
    const double offset_x = region_x + pad + (avail_w - span_x * scale_x) / 2.0 - min_x * scale_x;
    const double offset_y = region_y + pad + (avail_h - span_y * scale_y) / 2.0 - min_y * scale_y;

    std::vector<Point> screen(n);
    for (std::size_t i = 0; i < n; ++i) {
        screen[i] = Point{offset_x + layout[i].x * scale_x, offset_y + layout[i].y * scale_y};
    }

    double max_flow = 1e-9;
    for (std::size_t k = 0; k < nl; ++k) {
        max_flow = std::max(max_flow, std::max(std::abs(result.sf[k]), std::abs(result.st[k])));
    }

    const bool show_labels = n <= 70;
    const bool show_voltage = n <= 45;
    const double node_scale = n > 90 ? 0.55 : (n > 45 ? 0.75 : 1.0);

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\""
        << height << "\" viewBox=\"0 0 " << width << ' ' << height << "\" role=\"img\">\n";

    svg << "<defs>\n";
    svg << "  <linearGradient id=\"voltageRamp\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"0\">\n";
    for (int step = 0; step <= 8; ++step) {
        const double t = static_cast<double>(step) / 8.0;
        const double vm = 0.90 + 0.20 * t;
        svg << "    <stop offset=\"" << num(t * 100.0, 1) << "%\" stop-color=\""
            << hex(voltage_colour(vm)) << "\"/>\n";
    }
    svg << "  </linearGradient>\n";
    svg << "  <radialGradient id=\"backdrop\" cx=\"50%\" cy=\"38%\" r=\"78%\">\n";
    svg << "    <stop offset=\"0%\" stop-color=\"#15203a\"/>\n";
    svg << "    <stop offset=\"100%\" stop-color=\"#0b1220\"/>\n";
    svg << "  </radialGradient>\n";
    svg << "</defs>\n";

    svg << "<rect width=\"100%\" height=\"100%\" fill=\"url(#backdrop)\"/>\n";

    const std::string title =
        options.title.empty() ? net.data().name + "   |   " + std::to_string(n) + " buses, " +
                                    std::to_string(nl) + " branches"
                              : options.title;
    const std::string subtitle =
        options.subtitle.empty()
            ? std::string(to_string(result.method)) + "   |   " +
                  std::to_string(result.status.iterations) + " iterations   |   losses " +
                  num(result.p_loss_mw, 2) + " MW / " + num(result.q_loss_mvar, 2) + " MVAr"
            : options.subtitle;

    svg << "<text x=\"" << num(side_margin) << "\" y=\"34\" fill=\"#f1f5f9\" font-family=\""
        << kSans << "\" font-size=\"21\" font-weight=\"600\">" << title << "</text>\n";
    svg << "<text x=\"" << num(side_margin) << "\" y=\"56\" fill=\"#94a3b8\" font-family=\""
        << kMono << "\" font-size=\"12.5\">" << subtitle << "</text>\n";

    // ---- branches
    svg << "<g stroke-linecap=\"round\">\n";
    for (std::size_t k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        const Branch& branch = net.branch(k);
        const Point& a = screen[net.index_of(branch.from)];
        const Point& b = screen[net.index_of(branch.to)];
        const double flow = std::max(std::abs(result.sf[k]), std::abs(result.st[k]));
        const double relative = flow / max_flow;
        const bool rated = branch.rate_a > 0.0;
        const std::string colour = hex(loading_colour(result.loading[k], rated));
        const double thickness = 1.1 + 4.0 * std::sqrt(relative);

        svg << "  <line x1=\"" << num(a.x) << "\" y1=\"" << num(a.y) << "\" x2=\"" << num(b.x)
            << "\" y2=\"" << num(b.y) << "\" stroke=\"" << colour << "\" stroke-width=\""
            << num(thickness) << "\" opacity=\"0.8\"/>\n";
    }
    svg << "</g>\n";

    // ---- flow direction markers
    svg << "<g>\n";
    for (std::size_t k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        const double magnitude = std::abs(result.sf[k]);
        if (magnitude / max_flow < 0.06) {
            continue;
        }
        const Branch& branch = net.branch(k);
        Point a = screen[net.index_of(branch.from)];
        Point b = screen[net.index_of(branch.to)];
        if (result.sf[k].real() < 0.0) {
            std::swap(a, b);
        }
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double angle = std::atan2(dy, dx) * kRadToDeg;
        const double ax = a.x + dx * 0.58;
        const double ay = a.y + dy * 0.58;
        const bool rated = branch.rate_a > 0.0;
        const std::string colour = hex(loading_colour(result.loading[k], rated));
        svg << "  <polygon points=\"0,-3.4 8,0 0,3.4\" fill=\"" << colour
            << "\" opacity=\"0.95\" transform=\"translate(" << num(ax) << ',' << num(ay)
            << ") rotate(" << num(angle, 1) << ")\"/>\n";
    }
    svg << "</g>\n";

    // ---- buses
    svg << "<g font-family=\"" << kMono << "\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        const BusType type = net.bus_type(i);
        const bool connected = result.connected[i];
        const Rgb fill = voltage_colour(result.vm[i]);
        const double radius = (type == BusType::Slack ? 16.0 : (type == BusType::PV ? 13.0 : 10.5)) *
                              node_scale;

        const std::string ring = type == BusType::Slack ? "#f8fafc"
                                 : (type == BusType::PV ? "#cbd5e1" : "#1e293b");
        const double ring_width = type == BusType::Slack ? 2.6 : (type == BusType::PV ? 2.0 : 1.4);

        svg << "  <circle cx=\"" << num(screen[i].x) << "\" cy=\"" << num(screen[i].y)
            << "\" r=\"" << num(radius) << "\" fill=\"" << hex(fill) << "\" stroke=\"" << ring
            << "\" stroke-width=\"" << num(ring_width) << "\""
            << (connected ? "" : " opacity=\"0.35\"") << "/>\n";

        if (show_labels) {
            const std::string text_colour = luminance(fill) > 0.58 ? "#0b1220" : "#f8fafc";
            svg << "  <text x=\"" << num(screen[i].x) << "\" y=\"" << num(screen[i].y + 3.6)
                << "\" fill=\"" << text_colour << "\" font-size=\""
                << num(10.5 * node_scale + 1.0, 1) << "\" font-weight=\"600\" text-anchor=\"middle\">"
                << net.bus(i).id << "</text>\n";
        }
        if (show_voltage) {
            svg << "  <text x=\"" << num(screen[i].x) << "\" y=\""
                << num(screen[i].y + radius + 11.5)
                << "\" fill=\"#cbd5e1\" stroke=\"" << kHalo
                << "\" stroke-width=\"2.6\" stroke-linejoin=\"round\" paint-order=\"stroke\""
                << " font-size=\"9\" text-anchor=\"middle\">"
                << num(result.vm[i], 3) << "</text>\n";
        }
    }
    svg << "</g>\n";

    // ---- footer: legends and key figures
    const double footer_y = static_cast<double>(height) - footer_height;
    svg << "<line x1=\"" << num(side_margin) << "\" y1=\"" << num(footer_y + 2) << "\" x2=\""
        << num(width - side_margin) << "\" y2=\"" << num(footer_y + 2)
        << "\" stroke=\"#1e293b\" stroke-width=\"1\"/>\n";

    svg << "<g font-family=\"" << kSans << "\">\n";
    svg << "  <text x=\"" << num(side_margin) << "\" y=\"" << num(footer_y + 26)
        << "\" fill=\"#cbd5e1\" font-size=\"11.5\" font-weight=\"600\">bus voltage (pu)</text>\n";
    svg << "  <rect x=\"" << num(side_margin) << "\" y=\"" << num(footer_y + 34)
        << "\" width=\"232\" height=\"11\" rx=\"3\" fill=\"url(#voltageRamp)\"/>\n";
    for (int step = 0; step <= 4; ++step) {
        const double t = static_cast<double>(step) / 4.0;
        const double vm = 0.90 + 0.20 * t;
        svg << "  <text x=\"" << num(side_margin + 232.0 * t) << "\" y=\"" << num(footer_y + 62)
            << "\" fill=\"#94a3b8\" font-family=\"" << kMono
            << "\" font-size=\"10.5\" text-anchor=\"middle\">" << num(vm, 2) << "</text>\n";
    }

    svg << "  <text x=\"" << num(side_margin + 300.0) << "\" y=\"" << num(footer_y + 26)
        << "\" fill=\"#cbd5e1\" font-size=\"11.5\" font-weight=\"600\">branch loading</text>\n";
    const double samples[3] = {0.25, 0.70, 1.05};
    for (int s = 0; s < 3; ++s) {
        const double y = footer_y + 40.0 + 18.0 * static_cast<double>(s);
        svg << "  <line x1=\"" << num(side_margin + 300.0) << "\" y1=\"" << num(y) << "\" x2=\""
            << num(side_margin + 356.0) << "\" y2=\"" << num(y) << "\" stroke=\""
            << hex(loading_colour(samples[s], true)) << "\" stroke-width=\""
            << num(1.2 + 4.0 * std::sqrt(samples[s] / 1.2)) << "\" stroke-linecap=\"round\"/>\n";
        svg << "  <text x=\"" << num(side_margin + 366.0) << "\" y=\"" << num(y + 4)
            << "\" fill=\"#94a3b8\" font-family=\"" << kMono
            << "\" font-size=\"10.5\">" << num(samples[s] * 100.0, 0) << " %</text>\n";
    }
    svg << "  <text x=\"" << num(side_margin + 300.0) << "\" y=\"" << num(footer_y + 100)
        << "\" fill=\"#64748b\" font-size=\"10.5\" font-family=\"" << kMono
        << "\">line width and arrow size follow the power flow</text>\n";

    std::ostringstream stats;
    stats << "losses     " << std::setw(9) << num(result.p_loss_mw, 3) << " MW   " << std::setw(9)
          << num(result.q_loss_mvar, 3) << " MVAr\n";
    stats << "voltage    " << std::setw(9) << num(result.min_voltage(), 4) << "    " << std::setw(9)
          << num(result.max_voltage(), 4) << "    pu\n";
    stats << "mismatch   " << std::setw(9) << scientific(result.status.final_mismatch)
          << " pu   after " << result.status.iterations << " iterations";

    std::size_t line = 0;
    std::istringstream lines(stats.str());
    std::string text;
    while (std::getline(lines, text)) {
        svg << "  <text x=\"" << num(width - side_margin) << "\" y=\""
            << num(footer_y + 30.0 + 19.0 * static_cast<double>(line++))
            << "\" fill=\"#cbd5e1\" font-family=\"" << kMono
            << "\" font-size=\"11\" text-anchor=\"end\" xml:space=\"preserve\">" << text
            << "</text>\n";
    }
    svg << "</g>\n";

    svg << "</svg>\n";
    return svg.str();
}

}  // namespace pf
