"""Draw the diagrams for zhihu-answer.md from the plane equations and pivot maps."""

from itertools import combinations, permutations, product
from math import atan2, hypot
from pathlib import Path

from draw_figures import Figure, face_cells

OUT = Path(__file__).with_name("figures")
INK, MUTED, PALE = "#233847", "#516878", "#c5cfd5"
TEAL, OCHRE = "#007b80", "#b66014"
FACE = ("#e4eaf0", "#f1f4f6", "#fafbfc")
VERTICES = tuple(product((-.5, .5), repeat=3))
EDGES = [(a, b) for a, b in combinations(VERTICES, 2) if sum(x != y for x, y in zip(a, b)) == 1]


def project(p, cx, cy, scale):
    x, y, z = p
    return cx + scale * (.91*x - .68*y), cy + scale * (.32*x + .53*y - z)


def face_vertices(axis):
    points = []
    for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)):
        p = [.5] * 3
        p[(axis + 1) % 3], p[(axis + 2) % 3] = u, v
        points.append(tuple(p))
    return points


def hidden(a, b):
    return all(x == -.5 for x, y in zip(a, b) if x == y)


def header(fig, x, y, title):
    fig.label(x, y, title, 23, bold=True, color=INK)


def area(poly, axis):
    i, j = (axis+1) % 3, (axis+2) % 3
    return abs(sum(a[i]*b[j] - a[j]*b[i] for a, b in zip(poly, poly[1:]+poly[:1]))) / 2


def candidate():
    fig = Figure(656, 366)
    header(fig, 19, 36, "四刀在立方体表面的交线")
    header(fig, 350, 36, "右侧面  x = ½")
    draw = lambda p: project(p, 162, 171, 132)
    for axis in range(3):
        cells = list(face_cells(axis))
        assert abs(sum(area(poly, axis) for _, poly in cells)-1) < 1e-12
        for mask, poly in cells:
            fill = ("#76c4be" if mask == 3 else "#f6bf73" if mask == 1 else FACE[axis])
            fig.polygon([draw(p) for p in poly], fill, "#718791", 1.15)
        fig.polygon([draw(p) for p in face_vertices(axis)], "none", INK, 2)
    fig.line((286, 176), (327, 176), MUTED, 1.8)
    fig.polygon([(327, 171), (337, 176), (327, 181)], MUTED, MUTED, .8)
    left, top, width = 374, 73, 199
    def flat(p):
        return left + width*(.5-p[1]), top + width*(.5-p[2])
    for mask, poly in face_cells(0):
        fill = "#76c4be" if mask == 3 else "#f6bf73" if mask == 1 else "#f0f3f5"
        fig.polygon([flat(p) for p in poly], fill, "#536b77", 1.5)
    fig.polygon([flat(p) for p in face_vertices(0)], "none", INK, 2)
    # Each line is named inside the corner it removes from the highlighted block.
    fig.label(388, 95, "L₁", 16, color=INK)
    fig.label(545, 95, "L₃", 16, color=INK)
    fig.label(388, 259, "L₂", 16, color=INK)
    fig.label(545, 259, "L₀", 16, color=INK)
    fig.label(473, 178, "L₀、L₁ 为正", 17, "middle", color=INK, bold=True)
    fig.line((19, 307), (636, 307), "#d7dfe2", 1)
    for x, fill, label in [(24, "#76c4be", "一个大块的底面"), (328, "#f6bf73", "同一小块的三个底面三角形")]:
        fig.polygon([(x,328),(x+15,328),(x+15,343),(x,343)], fill, INK, .9)
        fig.label(x+25, 342, label, 18, color=INK)
    fig.save("answer-candidate")


def section(normal, level):
    def value(p):
        return sum(a*b for a,b in zip(p,normal))-level
    selected = [v for v in VERTICES if value(v)>0]
    points = []
    for a,b in EDGES:
        va, vb = value(a), value(b)
        if va*vb < 0:
            t = va/(va-vb)
            points.append(tuple(x+t*(y-x) for x,y in zip(a,b)))
    return selected, points


def vertex_types():
    fig = Figure(656, 494)
    for index,(name, count, n, h, size) in enumerate([
        ("一条棱",12,(1,1,0),.55,2), ("同面三点",24,(2,1,1),.45,3),
        ("一整个面",3,(1,0,0),.12,4), ("一个角星",4,(1,1,1),.12,4),
    ]):
        x,y = 18+(index%2)*329, 16+(index//2)*219
        header(fig, x, y+22, name)
        fig.label(x+278,y+22,f"{count} 种",18,"end",color=MUTED)
        draw = lambda p: project(p,x+139,y+123,83)
        selected,cut = section(n,h)
        assert len(selected) == size
        for a,b in EDGES:
            if hidden(a,b):
                fig.line(draw(a),draw(b), "#91a0a9", 1.3, "4 4")
        poly = [draw(p) for p in cut]
        cx,cy = (sum(q[k] for q in poly)/len(poly) for k in (0,1))
        poly.sort(key=lambda q: atan2(q[1]-cy,q[0]-cx))
        fig.polygon(poly, "#f9d69c", OCHRE, 1.45)
        for a,b in EDGES:
            if not hidden(a,b):
                fig.line(draw(a),draw(b), INK, 1.65)
        for a,b in EDGES:
            if a in selected and b in selected:
                fig.line(draw(a),draw(b), TEAL, 3)
        for v in VERTICES:
            active = v in selected
            fig.circle(*draw(v), 4.7 if active else 2.8, TEAL if active else "#fff",
                       "#064f53" if active else "#8496a1", 1.1)
    fig.line((18,224),(636,224),"#d7dfe2",1)
    fig.line((18,442),(636,442),"#d7dfe2",1)
    fig.circle(29,468,5.5,TEAL,TEAL,1)
    fig.label(44,475,"选中一侧的顶点",18,color=INK)
    fig.polygon([(320,461),(337,461),(337,476),(320,476)],"#f9d69c",OCHRE,1)
    fig.label(348,475,"截平面与立方体的交面",18,color=INK)
    fig.save("answer-vertex-types")


def arrow(fig, start, end, double=False):
    dx,dy = end[0]-start[0],end[1]-start[1]
    length = hypot(dx,dy)
    u,v = dx/length,dy/length
    a=(start[0]+14*u,start[1]+14*v)
    b=(end[0]-14*u,end[1]-14*v)
    fig.line(a,b,TEAL,2.4)
    def head(tip,i,j):
        x,y=tip
        fig.polygon([tip,(x-8*i-4*j,y-8*j+4*i),(x-8*i+4*j,y-8*j-4*i)],TEAL,TEAL,.5)
    head(b,u,v)
    if double:
        head(a,-u,-v)


def orbits():
    fig = Figure(656, 466)
    layouts = {
        "1032": ((-40,-20),(40,-20),(-40,22),(40,22)),
        "1230": ((-35,-20),(35,-20),(35,24),(-35,24)),
        "1200": ((-10,22),(35,-21),(55,24),(-70,-18)),
        "1000": ((0,0),(65,0),(-53,-24),(-53,25)),
        "1001": ((-22,0),(22,0),(-72,-23),(72,24)),
        "1002": ((-22,0),(40,0),(-60,29),(-75,-24)),
    }
    labels = ["两对互选","四点成环","三点环加一条入边","两个点同指一端","两个点分指两端","两步链指向一端"]
    representatives = set()
    for i, (word, layout) in enumerate(layouts.items()):
        out=list(map(int,word))
        orbit={tuple(perm[out[perm.index(j)]] for j in range(4)) for perm in permutations(range(4))}
        assert representatives.isdisjoint(orbit)
        representatives.update(orbit)
        x,y=111+215*(i%3), 112+190*(i//3)
        fig.label(x,y-63,word,24,"middle",color=INK,bold=True)
        points=[(x+a,y+b) for a,b in layout]
        for a,b in enumerate(out):
            pair=out[b]==a
            if not pair or a<b:
                arrow(fig,points[a],points[b],pair)
        for j,pt in enumerate(points):
            fig.circle(*pt,11,"#fff",INK,1.6)
            fig.label(pt[0],pt[1]+5.5,str(j),15,"middle",color=INK,bold=True)
        fig.label(x,y+61,labels[i],17,"middle",color=INK)
        fig.label(x,y+85,f"{len(orbit)} 种标号",16,"middle",color=MUTED)
    assert representatives == {f for f in product(range(4), repeat=4) if all(i != f[i] for i in range(4))}
    fig.line((18,211),(637,211),"#d7dfe2",1)
    fig.line((18,410),(637,410),"#d7dfe2",1)
    fig.label(24,443,"i → j 表示第 i 行选中第 j 列。双箭头表示互选。",18,color=INK)
    fig.save("answer-pivot-orbits")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    candidate()
    vertex_types()
    orbits()
