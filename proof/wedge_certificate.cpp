// Integer certificate for the following statement:
// If two central planes have unit normals with dot product 1/3,
// the part of the unit cube on which the two signs differ has volume <= 107/250.
// With --regular-tetrahedron, certify that four central planes with
// regular-tetrahedron normals have a block of volume at least 83/1000.
//
// Compile with: c++ -std=c++20 -O3 wedge_certificate.cpp -o wedge_certificate
#include <algorithm>
#include <array>
#include <charconv>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string_view>

using Integer = std::int64_t;
using WideInteger = __int128;

constexpr Integer grid = Integer{1} << 20;
constexpr int depth_limit = 13;
constexpr Integer volume_denominator = 24 * grid * grid;

struct Point {
    Integer x;
    Integer y;
    auto operator<=>(const Point &) const = default;
};

struct Line {
    Integer a;
    Integer b;
    Integer c;
};

Integer turn(Point a, Point b, Point c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

Integer grid_floor(WideInteger numerator, WideInteger denominator) {
    // The denominator is positive.
    numerator *= grid;
    return static_cast<Integer>(numerator / denominator - (numerator % denominator < 0));
}

// Work in [-1,1]^3 and on the part where u·x >= 0 and v·x <= 0.
// Return twice the area of an enclosing polygon in integer grid coordinates.
Integer twice_face_area(const std::array<Integer, 3> &u, const std::array<Integer, 3> &v, std::size_t axis,
                        int side) {
    const std::size_t j = (axis + 1) % 3;
    const std::size_t k = (axis + 2) % 3;
    const Line lines[] = {
        {1, 0, 1},
        {-1, 0, 1},
        {0, 1, 1},
        {0, -1, 1},
        {-u[j], -u[k], side * u[axis]},
        {v[j], v[k], -side * v[axis]},
    };
    Point points[60];
    int count = 0;
    for (int i = 0; i < 6; ++i) {
        for (int l = i + 1; l < 6; ++l) {
            const Line p = lines[i];
            const Line q = lines[l];
            WideInteger denominator = WideInteger(p.a) * q.b - WideInteger(p.b) * q.a;
            if (denominator == 0)
                continue;
            WideInteger x = WideInteger(p.c) * q.b - WideInteger(p.b) * q.c;
            WideInteger y = WideInteger(p.a) * q.c - WideInteger(p.c) * q.a;
            if (denominator < 0) {
                denominator = -denominator;
                x = -x;
                y = -y;
            }
            bool feasible = true;
            for (const Line boundary : lines) {
                if (boundary.a * x + boundary.b * y > boundary.c * denominator) {
                    feasible = false;
                    break;
                }
            }
            if (!feasible)
                continue;

            // Enclose this exact rational vertex in its grid square.
            const Integer left = grid_floor(x, denominator);
            const Integer bottom = grid_floor(y, denominator);
            for (Integer dx : {0, 1}) {
                for (Integer dy : {0, 1})
                    points[count++] = {left + dx, bottom + dy};
            }
        }
    }
    if (count == 0)
        return 0;

    std::sort(points, points + count);
    count = int(std::unique(points, points + count) - points);
    Point hull[120];
    int length = 0;
    for (int i = 0; i < count; ++i) {
        while (length >= 2 && turn(hull[length - 2], hull[length - 1], points[i]) <= 0) {
            --length;
        }
        hull[length++] = points[i];
    }
    const int lower_length = length;
    for (int i = count - 2; i >= 0; --i) {
        while (length > lower_length && turn(hull[length - 2], hull[length - 1], points[i]) <= 0) {
            --length;
        }
        hull[length++] = points[i];
    }

    Integer area = 0;
    for (int i = 1; i < length; ++i) {
        area += hull[i - 1].x * hull[i].y - hull[i - 1].y * hull[i].x;
    }
    return area;
}

using Vector = std::array<Integer, 3>;

Integer disagreement_numerator(const Vector &u, const Vector &v) {
    Integer area = 0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        for (int side : {-1, 1})
            area += twice_face_area(u, v, axis, side);
    }
    // Each face gives a pyramid of height 1. Double for the opposite sign
    // region and divide by the volume 8 of the working cube.
    return area;
}

// Use the quaternion (d,x,y,z) up to normalization. The columns below are
// multiplied by d*d+x*x+y*y+z*z.
std::array<Vector, 3> rotation_columns(Integer x, Integer y, Integer z, Integer d) {
    const std::array<Integer, 3> b = {
        d * d + x * x - y * y - z * z,
        2 * (x * y + d * z),
        2 * (x * z - d * y),
    };
    const std::array<Integer, 3> c = {
        2 * (x * y - d * z),
        d * d - x * x + y * y - z * z,
        2 * (y * z + d * x),
    };
    const std::array<Integer, 3> e = {
        2 * (x * z + d * y),
        2 * (y * z - d * x),
        d * d - x * x - y * y + z * z,
    };
    return {b, c, e};
}

Integer acute_pair_numerator(const std::array<Vector, 3> &columns, Integer cosine_weight,
                             Integer sine_weight) {
    const Vector &b = columns[0];
    const Vector &c = columns[1];
    Vector u, v;
    for (std::size_t i = 0; i < 3; ++i) {
        u[i] = cosine_weight * b[i] + sine_weight * c[i];
        v[i] = cosine_weight * b[i] - sine_weight * c[i];
    }
    return disagreement_numerator(u, v);
}

Integer regular_tetrahedron_numerator(const std::array<Vector, 3> &columns) {
    constexpr std::array<Integer, 4> sx = {1, 1, -1, -1};
    constexpr std::array<Integer, 4> sy = {1, -1, 1, -1};
    std::array<Vector, 4> normals;
    for (std::size_t i = 0; i < normals.size(); ++i) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            normals[i][axis] =
                sx[i] * columns[0][axis] + sy[i] * columns[1][axis] + sx[i] * sy[i] * columns[2][axis];
        }
    }

    // The sum of the six D(n_i, -n_j) is 3 minus the volume of the
    // union of the six 2+2 cells. At each point the count is either 2 or 3.
    Integer sum = 0;
    for (std::size_t i = 0; i < normals.size(); ++i) {
        for (std::size_t j = i + 1; j < normals.size(); ++j) {
            const Vector v = {-normals[j][0], -normals[j][1], -normals[j][2]};
            sum += disagreement_numerator(normals[i], v);
        }
    }
    return sum;
}

enum class Certificate { AcutePair, RationalAcutePair, SharpRationalAcutePair, RegularTetrahedron };

struct PairAnchor {
    Integer cosine = 100;
    Integer sine = 61;
    Integer upper_thousandths = 384;
};

struct Counts {
    Integer certified = 0;
    Integer excluded = 0;
    Integer subdivisions = 0;
    int deepest = 0;
};

// The cell's intervals are [(x-1)/d,(x+1)/d] and similarly for y and z.
// At depth <= 13, d <= 2^14 and all line coefficients have magnitude < 2^38.
// Intersection predicates and grid scaling then fit in signed 128-bit integers.
bool verify(Integer x, Integer y, Integer z, Integer d, int depth, Certificate certificate,
            const PairAnchor &anchor, Counts &counts) {
    const Integer ax = std::max<Integer>(0, std::abs(x) - 1);
    const Integer ay = std::max<Integer>(0, std::abs(y) - 1);
    const Integer az = std::max<Integer>(0, std::abs(z) - 1);
    if (std::max({ax, ay, az}) * 1000 > d * 415 || ax + ay + az > d) {
        ++counts.excluded;
        return true;
    }

    const auto columns = rotation_columns(x, y, z, d);
    bool certified;
    if (certificate == Certificate::AcutePair) {
        const Integer numerator = acute_pair_numerator(columns, 239, 169);
        // Rotation and rational approximation errors contribute 6/d and 1/30000.
        certified = WideInteger(numerator) * d * 30000 + WideInteger(volume_denominator) * (180000 + d) <=
                    WideInteger(12840) * d * volume_denominator;
    } else if (certificate == Certificate::RationalAcutePair ||
               certificate == Certificate::SharpRationalAcutePair) {
        const Integer numerator = acute_pair_numerator(columns, anchor.cosine, anchor.sine);
        // The sharp mode uses the verified .51 upper angular derivative.
        // Each normal turns by at most 2*sqrt(3)/d across this box, so the
        // two disagreements sum to less than 3.6/d. The legacy mode uses 6/d.
        const Integer rotation_error = certificate == Certificate::SharpRationalAcutePair ? 3600 : 6000;
        certified = WideInteger(numerator) * d * 1000 + WideInteger(volume_denominator) * rotation_error <=
                    WideInteger(anchor.upper_thousandths) * d * volume_denominator;
    } else {
        const Integer numerator = regular_tetrahedron_numerator(columns);
        // The count can change by at most 1 on the union of the four
        // sign-change regions. This gives a joint rotation error of 12/d.
        // Check numerator/volume_denominator + 12/d <= 1251/500.
        certified = WideInteger(numerator) * d * 500 + WideInteger(volume_denominator) * 6000 <=
                    WideInteger(1251) * d * volume_denominator;
    }
    if (certified) {
        ++counts.certified;
        counts.deepest = std::max(counts.deepest, depth);
        return true;
    }
    if (depth == depth_limit) {
        std::cerr << "Unresolved cell (x,y,z,d): " << x << ' ' << y << ' ' << z << ' ' << d << '\n';
        return false;
    }

    ++counts.subdivisions;
    for (int dx : {-1, 1}) {
        for (int dy : {-1, 1}) {
            for (int dz : {-1, 1}) {
                if (!verify(2 * x + dx, 2 * y + dy, 2 * z + dz, 2 * d, depth + 1, certificate, anchor,
                            counts)) {
                    return false;
                }
            }
        }
    }
    return true;
}

int main(int argc, char **argv) {
    Certificate certificate = Certificate::AcutePair;
    PairAnchor anchor;
    if (argc == 2 && std::string_view(argv[1]) == "--regular-tetrahedron") {
        certificate = Certificate::RegularTetrahedron;
    } else if (argc == 2 && std::string_view(argv[1]) == "--pair-61-100") {
        certificate = Certificate::RationalAcutePair;
    } else if (argc == 5 &&
               (std::string_view(argv[1]) == "--pair" || std::string_view(argv[1]) == "--sharp-pair" ||
                std::string_view(argv[1]) == "--sharp-pair-any-angle")) {
        const std::string_view mode(argv[1]);
        const auto read = [](std::string_view text, Integer &value) {
            auto result = std::from_chars(text.data(), text.data() + text.size(), value);
            return result.ec == std::errc{} && result.ptr == text.data() + text.size() && 0 < value &&
                   value <= 239;
        };
        std::string_view bound_text(argv[4]);
        auto bound = std::from_chars(bound_text.data(), bound_text.data() + bound_text.size(),
                                     anchor.upper_thousandths);
        if (!read(argv[2], anchor.cosine) || !read(argv[3], anchor.sine) ||
            (mode != "--sharp-pair-any-angle" && anchor.sine >= anchor.cosine) || bound.ec != std::errc{} ||
            bound.ptr != bound_text.data() + bound_text.size() || anchor.upper_thousandths <= 0 ||
            anchor.upper_thousandths >= 1000) {
            std::cerr << "Expected positive COSINE and SINE up to 239 and 0 < BOUND < 1000. "
                         "Use --sharp-pair-any-angle when SINE >= COSINE.\n";
            return 2;
        }
        certificate = mode == "--pair" ? Certificate::RationalAcutePair : Certificate::SharpRationalAcutePair;
    } else if (argc != 1) {
        std::cerr << "Usage: wedge_certificate [--regular-tetrahedron | --pair-61-100"
                     " | --pair COSINE SINE BOUND | --sharp-pair COSINE SINE BOUND"
                     " | --sharp-pair-any-angle COSINE SINE BOUND]\n";
        return 2;
    }
    Counts counts;
    // Left cube rotations permit max(|x/d|,|y/d|,|z/d|) <= sqrt(2)-1
    // and (|x|+|y|+|z|)/d <= 1. Conjugation by diagonal sign rotations
    // permits x,y >= 0. The two boxes below cover this part of the chart.
    const bool ok = verify(1, 1, -1, 4, 1, certificate, anchor, counts) &&
                    verify(1, 1, 1, 4, 1, certificate, anchor, counts);
    std::cout << (ok ? "certified" : "failed") << " leaves=" << counts.certified
              << " excluded=" << counts.excluded << " subdivided=" << counts.subdivisions
              << " depth=" << counts.deepest << '\n';
    return ok ? 0 : 1;
}
