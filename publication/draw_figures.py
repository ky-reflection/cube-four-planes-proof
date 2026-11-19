"""Generate the article's vector diagrams using only the Python standard library."""
from html import escape
from itertools import combinations, product
from math import sqrt
from pathlib import Path

OUT = Path(__file__).with_name("figures")
EDGE, INK, MUTED = "#27445a", "#18344a", "#526a7c"
TEAL, BLUE, GOLD = "#9ddad5", "#ccdeef", "#f9c777"
NORMALS = ((1, 1, 1), (sqrt(2), -1, -1), (-1, sqrt(2), -1), (-1, -1, sqrt(2)))
VERTICES = list(product((0, 1), repeat=3))


def number(x):
    return f"{x:.4f}"


class Figure:
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.parts = [f'<rect width="{width}" height="{height}" fill="#fff"/>']

    def line(self, a, b, color=EDGE, width=1.4, dash=None):
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(f'<path d="M{number(a[0])} {number(a[1])} L{number(b[0])} {number(b[1])}"'
                          f' fill="none" stroke="{color}" stroke-width="{width}"{dash_attr}/>')

    def polygon(self, points, color, stroke=EDGE, width=1.4):
        pts = " ".join(f"{number(x)},{number(y)}" for x, y in points)
        self.parts.append(f'<polygon points="{pts}" fill="{color}" stroke="{stroke}" '
                          f'stroke-width="{width}" stroke-linejoin="round"/>')

    def label(self, x, y, value, size=17, anchor="start", color=INK, bold=False):
        weight = 600 if bold else 400
        self.parts.append(f'<text x="{number(x)}" y="{number(y)}" font-size="{size}" '
                          f'font-weight="{weight}" fill="{color}" text-anchor="{anchor}" '
                          f'font-family="Hiragino Sans GB, Noto Sans CJK SC, sans-serif">{escape(value)}</text>')

    def circle(self, x, y, radius, color, stroke="white", width=1.6):
        self.parts.append(f'<circle cx="{number(x)}" cy="{number(y)}" r="{radius}" '
                          f'fill="{color}" stroke="{stroke}" stroke-width="{width}"/>')

    def save(self, name):
        self.parts.append("</svg>")
        tag = f'<svg xmlns="http://www.w3.org/2000/svg" width="{2*self.width}" height="{2*self.height}" viewBox="0 0 {self.width} {self.height}">'
        (OUT / (name + ".svg")).write_text(tag + "\n" + "\n".join(self.parts) + "\n")


def clip(poly, n, sign):
    """Clip a face polygon against sign * (n dot x) >= 0, interpolating on its edges."""
    if not poly:
        return []
    result = []
    for a, b in zip(poly, poly[1:] + poly[:1]):
        va, vb = (sign * sum(x * y for x, y in zip(n, q)) for q in (a, b))
        if va >= 0:
            result.append(a)
        if va * vb < 0:
            t = va / (va - vb)
            result.append(tuple(x + t * (y - x) for x, y in zip(a, b)))
    return result


def face_cells(axis):
    face = []
    for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)):
        q = [.5] * 3
        q[(axis + 1) % 3], q[(axis + 2) % 3] = u, v
        face.append(tuple(q))
    for mask in range(1, 15):
        poly = face
        for i, n in enumerate(NORMALS):
            poly = clip(poly, n, 1 if mask & (1 << i) else -1)
        if len(poly) < 3:
            continue
        p, q = (axis + 1) % 3, (axis + 2) % 3
        area2 = abs(sum(a[p]*b[q] - a[q]*b[p] for a, b in zip(poly, poly[1:] + poly[:1])))
        if area2 > 1e-10:
            yield mask, poly


def color(mask):
    if mask in (1, 14):
        return GOLD
    return TEAL if bin(mask).count("1") == 2 else BLUE


def isometric(p, x, y, scale):
    a, b, c = p
    return (x + scale * .8660254 * (a-b), y + scale * ((a+b)/2-c))


def candidate():
    f = Figure(620, 405)
    for face in (0, 1, 2):
        for mask, polygon in face_cells(face):
            f.polygon([isometric(p, 205, 202, 164) for p in polygon], color(mask), width=1.6)
    f.label(384, 112, "14 个空间块", 25, bold=True)
    for j, (shade, name, sub) in enumerate(((TEAL, "6 个两正块", "每块体积 m"),
                                              (BLUE, "6 个较大的单正或三正块", "每块体积 m"),
                                              (GOLD, "2 个较小的对顶块", "每块体积约 0.02145"))):
        y = 154 + 66*j
        f.circle(394, y-4, 7, shade, EDGE, .8)
        f.label(410, y, name, 15)
        f.label(410, y+22, sub, 14, color=MUTED)
    f.save("candidate")


def face():
    f = Figure(650, 387)
    def flat(p):
        return (183 + 275*p[1], 183 - 275*p[2])
    for mask, poly in face_cells(0):
        f.polygon([flat(p) for p in poly], color(mask), width=1.8)
    f.label(186, 193, "L₀、L₁ 为正", 16, "middle", bold=True)
    f.label(310, 70, "A", 13, "middle")
    f.label(60, 310, "B", 13, "middle")
    f.label(69, 69, "C", 13, "middle")
    f.label(305, 310, "D", 13, "middle")
    f.label(45, 359, "正方形边长为 1，y 轴向右，z 轴向上。", 15, color=MUTED)
    f.label(360, 89, "x = 1/2 面上的四个角", 20, bold=True)
    for i, (name, lengths) in enumerate((("A　小块", "1 − 1/√2，1 − 1/√2"),
                                         ("B", "1/2，1/2"), ("C、D", "1/2，1/√2"))):
        y = 132 + i*60
        f.label(362, y, name, 17, bold=True)
        f.label(362, y+23, "直角边　" + lengths, 15, color=MUTED)
    f.save("face-section")


def classes():
    f = Figure(660, 378)
    cases = [("棱", {(1,1,0),(1,1,1)}, 12),
             ("同一面上三个顶点", {(1,1,0),(1,1,1),(1,0,1)}, 24),
             ("整面", {p for p in VERTICES if p[0]}, 3),
             ("角星", {p for p in VERTICES if sum(p)>=2}, 4)]
    for i,(label,chosen,count) in enumerate(cases):
        x,y = 97+318*(i%2), 95+184*(i//2)
        def draw(p):
            a,b,c=(v-.5 for v in p)
            return (x+55*a-31*b, y+22*a+29*b-64*c)
        for a,b in combinations(VERTICES,2):
            if sum(s!=t for s,t in zip(a,b))==1:
                rear = not any(a) or not any(b)
                f.line(draw(a),draw(b),"#a6b4bd" if rear else EDGE,1.3,"4 3" if rear else None)
        for p in sorted(VERTICES,key=lambda v:sum(v)):
            px,py=draw(p)
            f.circle(px,py,5.3,"#008278" if p in chosen else "#ffffff", "#00665e" if p in chosen else "#8a9da9",1.4)
        f.label(x+69,y-3,label,17,bold=True)
        f.label(x+69,y+23,f"{count} 种无向分法",14,color=MUTED)
    f.save("vertex-types")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    candidate()
    face()
    classes()
