#define main pivot_cover_main
#include "affine_alternating_candidate_cover.cpp"
#undef main
int main() {
    int graph, pair;
    DirectionPairs cache;
    int last = -1;
    while (std::cin >> graph >> pair) {
        if (last != graph) {
            cache = DirectionPairs{};
            last = graph;
        }
        std::string s = std::to_string(graph);
        Pivots pivot{};
        Matrix c{};
        for (int i = 0; i < 4; ++i) {
            pivot[i] = s[i] - '0';
            for (int j = 0; j < 4; ++j)
                if (j != i) {
                    if (j == pivot[i])
                        c[i][j] = {grid, grid};
                    else
                        std::cin >> c[i][j].lo >> c[i][j].hi;
                }
        }
        auto q = cache.compute(c, pivot, pair);
        for (auto x : q.corner)
            std::cout << x.l << ' ' << x.h << ' ';
        for (int k = 0; k < 4; k++) {
            auto i = edges[pair][k / 2];
            int t = 0;
            Int radius = 0;
            for (int j = 0; j < 4; j++)
                if (j != i && j != pivot[i] && t++ == k % 2)
                    radius = (c[i][j].hi - c[i][j].lo) * (mass / grid) / 2;
            auto d = q.deriv[k / 2][k % 2], H = q.hessian[k / 2][k % 2];
            Int lip = ceil_ratio(Wide(d.h - d.l) * radius, 2 * mass);
            I errh = mul(divide(sq({radius, radius}), 2), H);
            I err = !q.has[k / 2][k % 2] ? I{-lip, lip}
                                         : I{std::max(-lip, -std::max<Int>(0, errh.h)),
                                             std::min(lip, -std::min<Int>(0, errh.l))};
            std::cout << err.l << ' ' << err.h << ' ';
        }
        std::cout << '\n';
    }
}
