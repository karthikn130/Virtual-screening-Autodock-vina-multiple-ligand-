# 🧬 Molecular Docking & Virtual Screening Suite

A standalone, integrated desktop application built in Python using **PySide6 (Qt)** to streamline molecular docking protocols, dynamic pocket visualization, and high-throughput virtual screening workflows. 

> ⚠️ **Platform Availability:** This application is compiled, optimized, and bundled **strictly for Linux environments** (tested heavily on Arch/Manjaro and Ubuntu distributions). 

---

## 🔬 Core Scientific Architecture

This software combines several open-source structural biology engines into a single graphical workbench:

* **AutoDock Vina Core:** Integrates native Python bindings for `vina` to execute rapid gradient-based empirical free energy force field calculations.
* **Open Babel Engine:** Handles automated on-the-fly 3D protonation states, structural file conversions (PDB/MOL2 to PDBQT), atom typing, and energy minimization protocols using the **MMFF94** force field.
* **PyMOL Visualizer Integration:** Connects seamlessly to render dynamic pocket selection, bounding box constraints, and real-time structural interaction vectors.

---

## 📂 Repository Directory Layout

```text
vina_compile/
├── main.py            # Application lifecycle entry point
├── gui.py            # PySide6 Qt GUI layout interfaces & windows
├── worker.py         # Multi-threaded QRunnable workers handling background calculations
└── visualiser.py     # Native rendering engine adapters (PyMOL interactions)
