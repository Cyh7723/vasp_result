import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
with (root/"cineb_energy_profile.csv").open() as handle:
    rows=list(csv.DictReader(handle))
x=[int(row["image"]) for row in rows]
y=[float(row["relative_to_initial_eV"]) for row in rows]
fig,ax=plt.subplots(figsize=(7.2,4.6),dpi=180)
ax.plot(x,y,color="#2457A6",marker="o",linewidth=1.8,markersize=6)
for index,energy in enumerate(y):
    ax.annotate(f"{energy:.3f}",(x[index],energy),textcoords="offset points",xytext=(0,8),ha="center",fontsize=9)
ax.set_xticks(x,[f"{i:02d}" for i in x])
ax.set_xlabel("CI-NEB image (00 = initial, 04 = final)")
ax.set_ylabel("Energy relative to initial (eV)")
ax.set_title("LiF V_Li(-) CI-NEB: 15-atom small-cell test")
ax.grid(axis="y",linestyle="--",alpha=0.35)
ax.spines[["top","right"]].set_visible(False)
ax.set_ylim(min(y)-0.1,max(y)+0.14)
fig.tight_layout()
fig.savefig(root/"cineb_energy_profile.png",bbox_inches="tight")
fig.savefig(root/"cineb_energy_profile.svg",bbox_inches="tight")
print("Saved energy profile PNG and SVG")
