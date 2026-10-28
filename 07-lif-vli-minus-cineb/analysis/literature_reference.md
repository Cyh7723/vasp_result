# 文献核对记录

来源：本地论文 `LiFNaF迁移能计算.pdf`（未随仓库分发）

论文：Handan Yildirim, Alper Kinaci, Maria K. Y. Chan, Jeffrey P. Greeley,
"First-Principles Analysis of Defect Thermodynamics and Ion Transport in
Inorganic SEI Compounds: LiF and NaF", ACS Applied Materials & Interfaces,
2015, 7, 18985–18996. DOI: 10.1021/acsami.5b02904.

已读取原文：

- 第3页（PDF页码）：PAW/PBE，LiF ENCUT=500 eV，5×5×5 Monkhorst–Pack，
  NEB 使用 64 原子超胞，自旋极化，最终力收敛 0.01 eV/Å。
- 同页 2.3 节：使用 3 个像；5 个像的测试给出相同过渡态几何和迁移势垒。
  方法正文未给出完整 INCAR，也未明确说明 LCLIMB 设置。
- 第7页表3：中性 Li 空位（V_Li^x）迁移势垒 0.73 eV；
  负电 Li 空位（V_Li'）迁移势垒 0.75 eV。
- 第8页正文再次给出中性 0.73 eV、负电 0.75 eV。

本次 15 原子测试不追求复现文献值。小超胞、简化 2×2×2 k 点和本次
VTST CI-NEB 优化方案均不同于完整文献生产计算；不得为接近 0.75 eV
而调整、挑选或修改输出。
