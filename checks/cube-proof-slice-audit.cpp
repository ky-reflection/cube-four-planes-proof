#include "affine_face_integral.hpp"
#include <iostream>
int main() {
    using namespace seven_integer;
    constexpr Bounds cube = {{{-space_grid / 2, space_grid / 2},
                              {-space_grid / 2, space_grid / 2},
                              {-space_grid / 2, space_grid / 2}}};
    std::array<Vector, 4> n{};
    std::array<Int, 4> h{};
    while (std::cin >> n[0][0]) {
        for (int i = 0; i < 4; i++) {
            for (int j = (i == 0); j < 3; j++)
                if (!(std::cin >> n[i][j]))
                    return 2;
            if (!(std::cin >> h[i]))
                return 2;
        }
        std::cout << volume<Bound::Lower>(n, h, cube) << ' ' << volume<Bound::Upper>(n, h, cube) << ' '
                  << cube_face_integral::volume<Bound::Lower>(n, h) << ' '
                  << cube_face_integral::volume<Bound::Upper>(n, h) << '\n';
    }
}
