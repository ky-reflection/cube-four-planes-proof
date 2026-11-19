"""Check local Taylor lower bounds in every parameter orthant using integer minors."""
from itertools import permutations
from pathlib import Path
import subprocess
import sys


def check_positive(matrix):
    """Sylvester criterion using fraction-free elimination without row exchanges."""
    last = 1
    for k in range(len(matrix)):
        pivot = matrix[k][k]
        assert pivot > 0, (k, pivot)
        for i in range(k + 1, len(matrix)):
            for j in range(i, len(matrix)):
                numerator = matrix[i][j] * pivot - matrix[i][k] * matrix[k][j]
                value, remainder = divmod(numerator, last)
                assert remainder == 0
                matrix[i][j] = matrix[j][i] = value
        last = pivot


def coefficients(binary, denominator, dimension, orthant=None):
    command = [str(Path(binary).resolve()), str(denominator)]
    if orthant is not None:
        command.append(str(orthant))
    run = subprocess.run(command,
                         capture_output=True, text=True, check=True)
    center = {}
    third = {}
    for line in run.stdout.splitlines():
        name, *columns = line.split()
        *indices, lo, hi = map(int, columns)
        assert lo <= hi
        if name == "H":
            assert len(indices)==2 and tuple(indices) not in center
            center[tuple(indices)] = (lo, hi)
        else:
            assert name == "T" and len(indices)==3 and tuple(indices) not in third
            for index in set(permutations(indices)):
                third[index] = (lo, hi)
    assert len(center)==dimension*(dimension+1)//2 and len(third)==dimension**3
    return center, third


def verify(binary, denominator, angular=False):
    dimension = 8 if angular else 12
    if not angular:
        center, third = coefficients(binary, denominator, dimension)
    factor=3*denominator
    for signs in range(1 << dimension):
        if angular:
            center, third = coefficients(binary, denominator, dimension, signs)
        negative=[(signs >> i) & 1 for i in range(dimension)]
        matrix=[[0]*dimension for _ in range(dimension)]
        for (i,j),(lo,hi) in center.items():
            opposite=negative[i]^negative[j]
            loss=0
            for k in range(dimension):
                l,u=third[i,j,k]
                loss+=min(0,-u if opposite^negative[k] else l)
            entry=factor*(hi if opposite else lo)+(-loss if opposite else loss)
            matrix[i][j]=matrix[j][i]=entry
        check_positive(matrix)
    print(f"radius 1/{denominator}: all {1 << dimension} Taylor lower forms are positive definite")


if __name__ == "__main__":
    if sys.argv[1] == "--angular":
        verify(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 35, angular=True)
    else:
        verify(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 100)
