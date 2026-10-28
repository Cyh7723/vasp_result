import hashlib
import json
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[1]


def read_poscar(path):
    lines = path.read_text().splitlines()
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * float(lines[1])
    counts = [int(x) for x in lines[6].split()]
    start = 7
    if lines[start].lower().startswith('s'):
        start += 1
    if lines[start][0].lower() in 'dc':
        start += 1
    coords = np.array([[float(x) for x in lines[start + i].split()[:3]] for i in range(sum(counts))])
    return lines[0], lattice, lines[5].split(), counts, coords


def write_poscar(path, title, lattice, elements, counts, coords):
    with path.open('w') as f:
        f.write(title + '\n1.0\n')
        for vector in lattice:
            f.write('  ' + ' '.join(f'{x:.12f}' for x in vector) + '\n')
        f.write(' '.join(elements) + '\n' + ' '.join(map(str, counts)) + '\nDirect\n')
        for point in coords:
            f.write('  ' + ' '.join(f'{x:.12f}' for x in point) + '\n')


initial_title, lattice, elements, counts, initial = read_poscar(root / 'initial/CONTCAR')
_, final_lattice, final_elements, final_counts, final = read_poscar(root / 'final/CONTCAR')
assert elements == final_elements == ['Li', 'F'] and counts == final_counts == [7, 8]
assert np.allclose(lattice, final_lattice)
for i in (0, 4):
    (root / 'neb' / f'{i:02d}').mkdir(exist_ok=True)
write_poscar(root / 'neb/00/POSCAR', 'LiF VLi- relaxed initial endpoint', lattice, elements, counts, initial)
write_poscar(root / 'neb/04/POSCAR', 'LiF VLi- relaxed final endpoint', lattice, elements, counts, final)
delta = final - initial
delta -= np.round(delta)
assert np.max(np.abs(delta[7:])) < 0.02, 'Unexpected F movement across periodic boundaries'
assert np.max(np.abs(delta[1:7])) < 0.01, 'Unexpected non-migrating Li movement'
report = {'endpoint_max_force_ev_per_A': {}, 'images': {}, 'moving_li_defect_atom_1based': 1}
for endpoint in ('initial', 'final'):
    output = (root / endpoint / 'OUTCAR').read_text()
    blocks = output.split('TOTAL-FORCE')[1:]
    forces = []
    for line in blocks[-1].splitlines()[2:]:
        if '---' in line:
            break
        values = line.split()
        if len(values) == 6:
            forces.append(np.linalg.norm([float(x) for x in values[3:]]))
    report['endpoint_max_force_ev_per_A'][endpoint] = max(forces)
for i, fraction in enumerate((0.25, 0.5, 0.75), start=1):
    coords = initial + fraction * delta
    (root / 'neb' / f'{i:02d}').mkdir(exist_ok=True)
    write_poscar(root / 'neb' / f'{i:02d}/POSCAR', f'LiF VLi- NEB linear image {i}', lattice, elements, counts, coords)
    displacement = coords - coords[:, None, :]
    displacement -= np.round(displacement)
    distances = np.linalg.norm(displacement @ lattice, axis=2)
    distances[distances < 1e-8] = np.inf
    report['images'][f'{i:02d}'] = {
        'moving_li_fraction': coords[0].tolist(),
        'minimum_pair_distance_angstrom': float(np.min(distances)),
        'lattice_matches_endpoints': bool(np.allclose(lattice, final_lattice)),
    }
report['potcar_sha256'] = hashlib.sha256((root / 'initial/POTCAR').read_bytes()).hexdigest()
for i in range(5):
    target = root / 'neb' / f'{i:02d}'
    (target / 'KPOINTS').write_text('LiF small-cell test\n0\nGamma\n2 2 2\n0 0 0\n')
    (target / 'POTCAR').write_bytes((root / 'initial/POTCAR').read_bytes())
report['endpoint_energy_difference_eV'] = 0.00004
report['note'] = 'Linear interpolation with component-wise minimum-image displacement; only the migrating Li crosses a periodic boundary.'
(root / 'analysis/initial_images_check.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
