import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[1]
neb = root / 'neb'
analysis = root / 'analysis'


def read_poscar(path):
    lines = path.read_text().splitlines()
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * float(lines[1])
    counts = [int(x) for x in lines[6].split()]
    start = 8 if lines[7][0].lower() in 'dc' else 7
    coords = np.array([[float(x) for x in lines[start + i].split()[:3]] for i in range(sum(counts))])
    return lines[0], lattice, lines[5].split(), counts, coords


def final_energy(path):
    values = re.findall(r'energy\(sigma->0\)\s*=\s*([-+\d.Ee]+)', path.read_text(errors='replace'))
    return float(values[-1]) if values else None


def min_distance(lattice, coords):
    delta = coords - coords[:, None, :]
    delta -= np.round(delta)
    distance = np.linalg.norm(delta @ lattice, axis=2)
    distance[distance < 1e-8] = np.inf
    return float(distance.min())


initial = final_energy(root / 'initial/OUTCAR')
final = final_energy(root / 'final/OUTCAR')
results = []
for i in range(5):
    directory = neb / f'{i:02d}'
    if i in (0, 4):
        energy = initial if i == 0 else final
        force = 0.0
        steps = 4
        text = (root / ('initial' if i == 0 else 'final') / 'OUTCAR').read_text()
        structure = directory / 'POSCAR'
    else:
        text = (directory / 'OUTCAR').read_text(errors='replace')
        energy = final_energy(directory / 'OUTCAR')
        force = [float(line.split()[-2]) for line in text.splitlines() if 'FORCES: max atom, RMS' in line][-1]
        steps = sum('F=' in line for line in (directory / 'OSZICAR').read_text().splitlines())
        structure = directory / 'CONTCAR'
    _, lattice, elements, counts, coords = read_poscar(structure)
    results.append({
        'image': f'{i:02d}', 'energy_sigma0_eV': energy,
        'relative_to_initial_eV': energy - initial,
        'max_neb_force_eV_per_A': force,
        'ionic_steps_this_stage': steps,
        'moving_li_fractional': coords[0].tolist(),
        'minimum_pair_distance_angstrom': min_distance(lattice, coords),
        'elements': elements, 'counts': counts,
        'electronic_ediff_stop_count': text.count('aborting loop because EDIFF is reached'),
        'normal_exit': 'General timing and accounting informations' in text,
        'lclimb_T': bool(re.search(r'NEB:\s+LCLIMB\s+T', text)),
        'ionic_converged': 'reached required accuracy' in text,
    })
maximum = max(item['energy_sigma0_eV'] for item in results)
barrier_forward = maximum - initial
barrier_reverse = maximum - final
with (analysis / 'cineb_energy_profile.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['image', 'energy_sigma0_eV', 'relative_to_initial_eV', 'max_neb_force_eV_per_A', 'moving_li_fractional', 'minimum_pair_distance_angstrom'])
    writer.writeheader()
    for row in results:
        writer.writerow({key: (json.dumps(value) if isinstance(value, list) else value) for key, value in row.items() if key in writer.fieldnames})
summary = {
    'status': 'converged_cineb_test',
    'model': 'LiF 2x2x2 fcc primitive-cell supercell; 7 Li + 8 F; negatively charged Li vacancy',
    'nelect': 78,
    'encut_eV': 500,
    'kpoints': 'Gamma-centered 2x2x2',
    'spin_polarized': True,
    'intermediate_images': 3,
    'mpi_ranks': 6,
    'ranks_per_image': 2,
    'omp_threads_per_rank': 1,
    'vasp_executable': 'User-provided VTST-enabled VASP 6.5.1',
    'vtst_version': '5.1, (04/25/25)',
    'lclimb_T_in_output': all(row['lclimb_T'] for row in results[1:4]),
    'final_electronic_ediff_stop_count': {row['image']: row['electronic_ediff_stop_count'] for row in results[1:4]},
    'all_three_images_ediff_each_step': all(row['electronic_ediff_stop_count'] == 7 for row in results[1:4]),
    'all_three_images_force_converged': all(row['ionic_converged'] and row['max_neb_force_eV_per_A'] <= 0.03 for row in results[1:4]),
    'max_final_neb_force_eV_per_A': max(row['max_neb_force_eV_per_A'] for row in results[1:4]),
    'cineb_steps': 7,
    'endpoint_max_force_ev_per_A': json.loads((analysis / 'initial_images_check.json').read_text())['endpoint_max_force_ev_per_A'],
    'endpoint_initial_energy_sigma0_eV': initial,
    'endpoint_final_energy_sigma0_eV': final,
    'endpoint_energy_difference_eV': final - initial,
    'forward_barrier_eV': barrier_forward,
    'reverse_barrier_eV': barrier_reverse,
    'highest_image': max(results, key=lambda row: row['energy_sigma0_eV'])['image'],
    'energy_definition': 'OUTCAR energy(sigma->0), endpoints from separate relaxations with same electronic parameters',
    'profile': results,
    'pre_relax': {'first_stage_steps': 5, 'first_stage_wall_seconds': 384.96,
                  'continued_stage_steps': 8, 'continued_stage_wall_seconds': 446.17,
                  'final_max_neb_force_eV_per_A': 0.047415},
    'cineb_wall_seconds': 332.44,
    'cineb_seconds_per_ionic_step': 332.44 / 7,
    'literature_reference': {'neutral_vacancy_barrier_eV': 0.73, 'negative_vacancy_barrier_eV': 0.75,
                              'caution': 'Small-cell educational test; not expected to reproduce paper values.'},
}
(analysis / 'cineb_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x = list(range(5))
    energies = [row['relative_to_initial_eV'] for row in results]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=180)
    ax.plot(x, energies, color='#2457A6', marker='o', linewidth=1.8, markersize=6)
    for index, energy in enumerate(energies):
        ax.annotate(f'{energy:.3f}', (x[index], energy), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=9)
    ax.set_xticks(x, [f'{i:02d}' for i in x])
    ax.set_xlabel('CI-NEB image (00 = initial, 04 = final)')
    ax.set_ylabel('Energy relative to initial (eV)')
    ax.set_title('LiF $V_{Li}^{-}$ CI-NEB: 15-atom small-cell test')
    ax.grid(axis='y', linestyle='--', alpha=0.35)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(min(energies) - 0.1, max(energies) + 0.14)
    fig.tight_layout()
    fig.savefig(analysis / 'cineb_energy_profile.png', bbox_inches='tight')
    plt.close(fig)
    print('Energy plot written')
except ImportError:
    print('matplotlib unavailable; CSV and JSON are complete')
