import os
import time
from PySide6.QtCore import QThread, Signal
from openbabel import pybel
from vina import Vina

class DockingWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    finished_signal = Signal(str)

    def __init__(self, receptor_path, ligand_input_path, grid_params):
        super().__init__()
        self.receptor_path = receptor_path
        self.ligand_input_path = ligand_input_path
        self.grid_params = grid_params
        self.is_running = True

    def run(self):
        try:
            output_dir = os.path.join(os.getcwd(), "docking_outputs")
            os.makedirs(output_dir, exist_ok=True)

            self.log_signal.emit("⏳ Preparing structural receptor pocket space coordinates...")
            rec_name = os.path.basename(self.receptor_path).split('.')[0]
            receptor_pdbqt = os.path.join(output_dir, f"{rec_name}_prepared.pdbqt")

            # Receptor Format Conversion
            if not self.receptor_path.lower().endswith('.pdbqt'):
                mol = next(pybel.readfile(self.receptor_path.split('.')[-1].lower(), self.receptor_path))
                mol.addh()
                mol.write("pdbqt", receptor_pdbqt, overwrite=True)
            else:
                receptor_pdbqt = self.receptor_path

            # Parse Ligand Inputs
            raw_ligand_files = []
            if os.path.isdir(self.ligand_input_path):
                for f in os.listdir(self.ligand_input_path):
                    if f.lower().endswith(('.sdf', '.mol2', '.pdbqt')):
                        raw_ligand_files.append(os.path.join(self.ligand_input_path, f))
            else:
                if self.ligand_input_path.lower().endswith(('.sdf', '.mol2', '.pdbqt')):
                    raw_ligand_files.append(self.ligand_input_path)

            total_ligands = len(raw_ligand_files)
            if total_ligands == 0:
                self.finished_signal.emit("❌ Stop: No valid ligand molecules detected.")
                return

            cores = self.grid_params.get('cores_per_test_system', os.cpu_count() or 4)
            self.log_signal.emit(f"🚀 Screening pipeline initialized for {total_ligands} compounds utilizing {cores} CPU cores...")

            for idx, lig_path in enumerate(raw_ligand_files):
                if not self.is_running:
                    break

                lig_name = os.path.basename(lig_path).split('.')[0]
                self.log_signal.emit(f"\n🧬 Processing Molecule [{idx+1}/{total_ligands}]: {lig_name}")
                
                lig_ext = lig_path.split('.')[-1].lower()
                target_pdbqt = os.path.join(output_dir, f"{lig_name}_prepared.pdbqt")

                # On-the-Fly Optimization & Conversion
                if lig_ext != "pdbqt":
                    mol = next(pybel.readfile(lig_ext, lig_path))
                    mol.addh()
                    mol.make3D(forcefield='mmff94', steps=50)
                    
                    ff = pybel._forcefields['mmff94']
                    ff.Setup(mol.OBMol)
                    initial_energy = ff.Energy()
                    
                    mol.localopt(forcefield='mmff94', steps=150)
                    minimized_energy = ff.Energy()
                    
                    self.log_signal.emit(f"   • Initial Energy  : {initial_energy:.4f} kJ/mol")
                    self.log_signal.emit(f"   • Minimized Energy: {minimized_energy:.4f} kJ/mol (ΔE: {minimized_energy - initial_energy:.4f})")
                    
                    mol.write("pdbqt", target_pdbqt, overwrite=True)
                    current_lig_pdbqt = target_pdbqt
                else:
                    current_lig_pdbqt = lig_path

                # Call the native Python Vina module bindings directly
                start_clock = time.time()
                try:
                    v = Vina(sf_name='vina', cpu=cores, verbosity=1)
                    v.set_receptor(receptor_pdbqt)
                    v.set_ligand_from_file(current_lig_pdbqt)
                    v.compute_vina_maps(center=self.grid_params['center'], box_size=self.grid_params['size'])
                    
                    self.log_signal.emit(f"   • Optimization grid locked. Computing global search...")
                    v.dock(exhaustiveness=self.grid_params.get('exhaustiveness', 8), 
                           n_poses=self.grid_params.get('n_poses', 5))

                    output_poses = os.path.join(output_dir, f"{lig_name}_results.pdbqt")
                    v.write_poses(output_poses, n_poses=self.grid_params.get('n_poses', 5), overwrite=True)
                    
                    elapsed = time.time() - start_clock
                    self.log_signal.emit(f"   • Success: Configuration completed in {elapsed:.2f}s.")
                except Exception as docking_err:
                    self.log_signal.emit(f"   • Engine Error: {str(docking_err)}")

                # Cleanup temporary raw conversion intermediates
                if lig_ext != "pdbqt" and os.path.exists(target_pdbqt):
                    try: os.remove(target_pdbqt)
                    except Exception: pass

                prog_pct = int(10 + (90 * (idx + 1) / total_ligands))
                self.progress_signal.emit(prog_pct)

            self.progress_signal.emit(100)
            self.finished_signal.emit("🎉 Complete! Screening pipeline finished cleanly.")
            
        except Exception as e:
            self.finished_signal.emit(f"❌ Execution breakdown: {str(e)}")

    def stop(self):
        self.is_running = False