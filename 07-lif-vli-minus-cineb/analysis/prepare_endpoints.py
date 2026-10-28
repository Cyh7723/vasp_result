import json
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[1]


def read_poscar(path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    start = 7
    selective = False
    if lines[start].lower().startswith('s'):
        selective = True
        start += 1
    if lines[start][0].lower() in 'dc':
        start += 1
    coords = np.array([[float(x) for x in lines[start + i].split()[:3]] for i in range(sum(counts))])
    return lines[0], lattice, elements, counts, coords, selective


def write_poscar(path, title, lattice, elements, counts, coords):
    with path.open('w') as f:
        f.write(title + '\n1.0\n')
        for vector in lattice:
            f.write('  ' + ' '.join(f'{x:.12f}' for x in vector) + '\n')
        f.write(' '.join(elements) + '\n')
        f.write(' '.join(str(x) for x in counts) + '\nDirect\n')
        for point in coords:
            f.write('  ' + ' '.join(f'{x:.12f}' for x in point) + '\n')


title, lattice, elements, counts, full, _ = read_poscar(root / 'reference/complete_16atom_POSCAR')
assert elements == ['Li', 'F'] and counts == [8, 8]
vacancy_index = 0
vacancy = full[vacancy_index].copy()
remaining_li = np.delete(full[:8], vacancy_index, axis=0)
delta = remaining_li - vacancy
delta -= np.round(delta)
distances = np.linalg.norm(delta @ lattice, axis=1)
nearest_original_li_index = int(np.argmin(distances))
moving_li_index = nearest_original_li_index + 1
start_li = full[moving_li_index].copy()
displacement = vacancy - start_li
displacement -= np.round(displacement)
assert np.linalg.norm(displacement @ lattice) == distances[nearest_original_li_index]
initial = np.vstack([np.delete(full[:8], vacancy_index, axis=0), full[8:]])
final = initial.copy()
final[nearest_original_li_index] = vacancy
assert np.allclose(final[nearest_original_li_index], [0, 0, 0])
write_poscar(root / 'initial/POSCAR', 'LiF VLi- initial: vacancy at original Li 1', lattice, ['Li', 'F'], [7, 8], initial)
write_poscar(root / 'final/POSCAR', 'LiF VLi- final: original Li 2 moved into vacancy', lattice, ['Li', 'F'], [7, 8], final)
write_poscar(root / 'initial/KPOINTS', 'LiF small-cell test', np.eye(3), [], [], np.empty((0, 3))) if False else None
(root / 'initial/KPOINTS').write_text('LiF small-cell test\n0\nGamma\n2 2 2\n0 0 0\n')
(root / 'final/KPOINTS').write_text('LiF small-cell test\n0\nGamma\n2 2 2\n0 0 0\n')
metadata = {
    'vacancy_original_full_atom_1based': vacancy_index + 1,
    'vacancy_fractional_coordinate': vacancy.tolist(),
    'moving_li_original_full_atom_1based': moving_li_index + 1,
    'moving_li_defect_atom_1based': nearest_original_li_index + 1,
    'start_fractional_coordinate': start_li.tolist(),
    'end_fractional_coordinate': vacancy.tolist(),
    'minimum_image_displacement_fractional': displacement.tolist(),
    'migration_distance_angstrom': float(np.linalg.norm(displacement @ lattice)),
    'candidate_neighbor_distances_angstrom': distances.tolist(),
    'atom_order': 'Delete original full Li 1; all other atom identities and ordering retained. Move original full Li 2 to vacancy for final.',
    'lattice_identical_for_all_endpoints_and_images': True,
    'nelect': 7 * 3 + 8 * 7 + 1,
    'potcar': [{'element': 'Li', 'TITEL': 'PAW_PBE Li_sv 10Sep2004', 'ZVAL': 3.0, 'ENMAX': 499.034},
               {'element': 'F', 'TITEL': 'PAW_PBE F 08Apr2002', 'ZVAL': 7.0, 'ENMAX': 400.0}],
    'vasp_executable': 'User-provided VTST-enabled VASP 6.5.1',
    'mpi_command': 'mpirun --oversubscribe -np 6 (3 images, 2 ranks per image)',
}
(root / 'analysis/endpoint_mapping.json').write_text(json.dumps(metadata, indent=2) + '\n')
print(json.dumps(metadata, indent=2))
