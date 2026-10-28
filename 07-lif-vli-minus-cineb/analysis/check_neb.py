import json
import re
import sys
from pathlib import Path

import numpy as np

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'neb'
result = {'root': str(root), 'images': {}}


def read_poscar(path):
    lines = path.read_text().splitlines()
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * float(lines[1])
    start = 8 if lines[7][0].lower() in 'dc' else 7
    counts = [int(x) for x in lines[6].split()]
    coords = np.array([[float(x) for x in lines[start + i].split()[:3]] for i in range(sum(counts))])
    return lattice, coords


for i in range(1, 4):
    directory = root / f'{i:02d}'
    outcar = directory / 'OUTCAR'
    if not outcar.exists():
        result['images'][f'{i:02d}'] = {'status': 'not_started'}
        continue
    text = outcar.read_text(errors='replace')
    energies = []
    for line in (directory / 'OSZICAR').read_text().splitlines():
        fields = line.split()
        if 'E0=' in fields:
            energies.append((int(fields[0]), float(fields[fields.index('E0=') + 1])))
    forces = [float(line.split()[-2]) for line in text.splitlines() if 'FORCES: max atom, RMS' in line]
    electronic_stop = text.count('aborting loop because EDIFF is reached')
    poscar = directory / 'CONTCAR' if (directory / 'CONTCAR').exists() else directory / 'POSCAR'
    lattice, coords = read_poscar(poscar)
    delta = coords - coords[:, None, :]
    delta -= np.round(delta)
    distances = np.linalg.norm(delta @ lattice, axis=2)
    distances[distances < 1e-8] = np.inf
    elapsed = re.findall(r'Elapsed time \(sec\):\s*([\d.]+)', text)
    result['images'][f'{i:02d}'] = {
        'ionic_steps': len(energies),
        'energies_E0_eV': energies,
        'latest_neb_max_force_eV_per_A': forces[-1] if forces else None,
        'max_neb_force_history_eV_per_A': forces,
        'electronic_ediff_stop_count': electronic_stop,
        'ionic_converged': 'reached required accuracy' in text,
        'normal_exit': 'General timing and accounting informations' in text,
        'vtst_initialized': 'VTST: version' in text,
        'lclimb_T': bool(re.search(r'NEB:\s+LCLIMB\s+T', text)),
        'lclimb_F': bool(re.search(r'NEB:\s+LCLIMB\s+F', text)),
        'moving_li_fractional': coords[0].tolist(),
        'minimum_pair_distance_angstrom': float(np.min(distances)),
        'elapsed_seconds': float(elapsed[-1]) if elapsed else None,
    }
result['max_latest_neb_force_eV_per_A'] = max((v['latest_neb_max_force_eV_per_A'] for v in result['images'].values() if v.get('latest_neb_max_force_eV_per_A') is not None), default=None)
print(json.dumps(result, indent=2))
(root.parent / 'analysis' / f'{root.name}_status.json').write_text(json.dumps(result, indent=2) + '\n')
