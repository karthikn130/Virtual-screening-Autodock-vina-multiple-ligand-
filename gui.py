import os
import csv
from PySide6.QtWidgets import (QMainWindow, QWidget, QTabWidget, QVBoxLayout, 
                              QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                              QFileDialog, QFormLayout, QProgressBar, QTextEdit,
                              QTableWidget, QTableWidgetItem, QHeaderView, QMessageBox)
from PySide6.QtCore import Qt
from worker import DockingWorker
from visualizer import launch_pymol

class DockingApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Custom Molecular Docking Suite")
        self.setGeometry(100, 100, 850, 600)
        
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        self.init_input_tab()
        self.init_grid_tab()
        self.init_execution_tab()
        self.init_results_tab()
        self.init_citation_tab()

    def init_input_tab(self):
        tab = QWidget()
        layout = QFormLayout()
        
        self.receptor_input = QLineEdit()
        btn_browse_rec = QPushButton("Browse")
        btn_browse_rec.clicked.connect(lambda: self.browse_file(self.receptor_input, "Protein Files (*.pdb *.cif *.pdbqt)"))
        row_rec = QHBoxLayout()
        row_rec.addWidget(self.receptor_input)
        row_rec.addWidget(btn_browse_rec)
        layout.addRow(QLabel("<b>Protein Receptor:</b>"), row_rec)
        
        self.ligand_input = QLineEdit()
        btn_browse_lig_file = QPushButton("File")
        btn_browse_lig_file.clicked.connect(lambda: self.browse_file(self.ligand_input, "Ligand Files (*.sdf *.mol2 *.pdbqt)"))
        btn_browse_lig_dir = QPushButton("Folder")
        btn_browse_lig_dir.clicked.connect(self.browse_directory)
        
        row_lig = QHBoxLayout()
        row_lig.addWidget(self.ligand_input)
        row_lig.addWidget(btn_browse_lig_file)
        row_lig.addWidget(btn_browse_lig_dir)
        layout.addRow(QLabel("<b>Ligand File/Folder:</b>"), row_lig)
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "1. Input Selection")

    def browse_file(self, line_edit, file_filter):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", file_filter)
        if file_path:
            line_edit.setText(file_path)

    def browse_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Ligand Folder")
        if dir_path:
            self.ligand_input.setText(dir_path)

    def init_grid_tab(self):
        tab = QWidget()
        layout = QFormLayout()
        layout.addRow(QLabel("<h3>Define Search Space (Grid Box)</h3>"))
        
        from PySide6.QtWidgets import QCheckBox
        
        self.chk_auto_center_full = QCheckBox("Auto-Center on Full Protein")
        self.chk_auto_center_full.clicked.connect(self.handle_center_toggles)
        layout.addRow(self.chk_auto_center_full)
        
        self.chk_auto_center_pocket = QCheckBox("Auto-Center on Specific Pocket (Residues)")
        self.chk_auto_center_pocket.clicked.connect(self.handle_center_toggles)
        layout.addRow(self.chk_auto_center_pocket)
        
        self.txt_residues = QLineEdit()
        self.txt_residues.setPlaceholderText("e.g., GLY143, CYS145, HIS163")
        self.txt_residues.setEnabled(False)
        self.txt_residues.editingFinished.connect(self.process_pocket_calculation)
        layout.addRow(QLabel("<b>Target Pocket Residues:</b>"), self.txt_residues)
        
        layout.addRow(QLabel("<hr style='border: 0; border-top: 1px dashed #ccc;'>"))
        
        self.chk_auto_size_full = QCheckBox("Auto-Size to Encompass Full Protein")
        self.chk_auto_size_full.clicked.connect(self.handle_size_toggles)
        layout.addRow(self.chk_auto_size_full)
        
        layout.addRow(QLabel("<br>"))

        self.cx = QLineEdit("0.0")
        self.cy = QLineEdit("0.0")
        self.cz = QLineEdit("0.0")
        grid_center = QHBoxLayout()
        grid_center.addWidget(QLabel("X:"))
        grid_center.addWidget(self.cx)
        grid_center.addWidget(QLabel("Y:"))
        grid_center.addWidget(self.cy)
        grid_center.addWidget(QLabel("Z:"))
        grid_center.addWidget(self.cz)
        layout.addRow(QLabel("<b>Grid Center:</b>"), grid_center)
        
        self.sx = QLineEdit("20.0")
        self.sy = QLineEdit("20.0")
        self.sz = QLineEdit("20.0")
        grid_size = QHBoxLayout()
        grid_size.addWidget(QLabel("Size X:"))
        grid_size.addWidget(self.sx)
        grid_size.addWidget(QLabel("Size Y:"))
        grid_size.addWidget(self.sy)
        grid_size.addWidget(QLabel("Size Z:"))
        grid_size.addWidget(self.sz)
        layout.addRow(QLabel("<b>Grid Size (Å):</b>"), grid_size)
        
        self.lbl_volume = QLabel("<b>Calculated Volume:</b> 8000.000 Å³")
        self.lbl_volume.setStyleSheet("color: #5cb85c; font-size: 13px; font-weight: bold; background-color: #f4f9f4; padding: 6px; border-radius: 4px;")
        layout.addRow(self.lbl_volume)
        
        self.sx.textChanged.connect(self.update_volume_display)
        self.sy.textChanged.connect(self.update_volume_display)
        self.sz.textChanged.connect(self.update_volume_display)
        
        layout.addRow(QLabel("<hr>"))
        layout.addRow(QLabel("<h3>Search Engine Settings</h3>"))
        
        self.txt_exhaustiveness = QLineEdit("8")
        self.txt_n_poses = QLineEdit("5")
        
        total_cores = os.cpu_count() or 4
        layout.addRow(QLabel(f"<b>Hardware Profiles:</b> Total System CPU Cores Detected: {total_cores}"))
        layout.addRow(QLabel("<b>Exhaustiveness:</b>"), self.txt_exhaustiveness)
        layout.addRow(QLabel("<b>Number of Poses to Save:</b>"), self.txt_n_poses)
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "2. Grid Parameters")

    def calculate_protein_center(self, file_path):
        x_coords, y_coords, z_coords = [], [], []
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM")):
                        x_coords.append(float(line[30:38].strip()))
                        y_coords.append(float(line[38:46].strip()))
                        z_coords.append(float(line[46:54].strip()))
            if not x_coords: return None
            return (sum(x_coords)/len(x_coords), sum(y_coords)/len(y_coords), sum(z_coords)/len(z_coords))
        except Exception: return None

    def calculate_protein_dimensions(self, file_path):
        x_coords, y_coords, z_coords = [], [], []
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM")):
                        x_coords.append(float(line[30:38].strip()))
                        y_coords.append(float(line[38:46].strip()))
                        z_coords.append(float(line[46:54].strip()))
            if not x_coords: return None
            return (max(x_coords)-min(x_coords)+10.0, max(y_coords)-min(y_coords)+10.0, max(z_coords)-min(z_coords)+10.0)
        except Exception: return None

    def handle_center_toggles(self):
        receptor_path = self.receptor_input.text()
        sender = self.sender()
        if sender.isChecked() and (not receptor_path or not os.path.exists(receptor_path)):
            QMessageBox.warning(self, "Input Missing", "Please select a valid Protein Receptor on Tab 1 first.")
            sender.setChecked(False)
            return
        if sender == self.chk_auto_center_full and self.chk_auto_center_full.isChecked():
            self.chk_auto_center_pocket.setChecked(False)
            self.txt_residues.setEnabled(False)
            center = self.calculate_protein_center(receptor_path)
            if center:
                self.cx.setText(f"{center[0]:.3f}"); self.cy.setText(f"{center[1]:.3f}"); self.cz.setText(f"{center[2]:.3f}")
                self.cx.setReadOnly(True); self.cy.setReadOnly(True); self.cz.setReadOnly(True)
        elif sender == self.chk_auto_center_pocket and self.chk_auto_center_pocket.isChecked():
            self.chk_auto_center_full.setChecked(False)
            self.txt_residues.setEnabled(True)
            self.cx.setReadOnly(True); self.cy.setReadOnly(True); self.cz.setReadOnly(True)
            self.process_pocket_calculation()
        if not self.chk_auto_center_full.isChecked() and not self.chk_auto_center_pocket.isChecked():
            self.txt_residues.setEnabled(False)
            self.cx.setReadOnly(False); self.cy.setReadOnly(False); self.cz.setReadOnly(False)

    def handle_size_toggles(self):
        receptor_path = self.receptor_input.text()
        if self.chk_auto_size_full.isChecked():
            if not receptor_path or not os.path.exists(receptor_path):
                QMessageBox.warning(self, "Input Missing", "Please select a valid Protein Receptor on Tab 1 first.")
                self.chk_auto_size_full.setChecked(False)
                return
            dims = self.calculate_protein_dimensions(receptor_path)
            if dims:
                self.sx.setText(f"{dims[0]:.3f}"); self.sy.setText(f"{dims[1]:.3f}"); self.sz.setText(f"{dims[2]:.3f}")
                self.sx.setReadOnly(True); self.sy.setReadOnly(True); self.sz.setReadOnly(True)
        else:
            self.sx.setReadOnly(False); self.sy.setReadOnly(False); self.sz.setReadOnly(False)
        self.update_volume_display()

    def update_volume_display(self):
        try:
            volume = float(self.sx.text()) * float(self.sy.text()) * float(self.sz.text())
            if volume > 27000:
                self.lbl_volume.setText(f"⚠️ Calculated Volume: {volume:,.3f} Å³ (Large box - consider high exhaustiveness)")
                self.lbl_volume.setStyleSheet("color: #d9534f; font-size: 13px; font-weight: bold; background-color: #fdf7f7; padding: 6px; border-radius: 4px;")
            else:
                self.lbl_volume.setText(f"✅ Calculated Volume: {volume:,.3f} Å³ (Optimal size)")
                self.lbl_volume.setStyleSheet("color: #5cb85c; font-size: 13px; font-weight: bold; background-color: #f4f9f4; padding: 6px; border-radius: 4px;")
        except ValueError:
            self.lbl_volume.setText("❌ Calculated Volume: Invalid coordinate entries.")

    def init_execution_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        self.btn_run = QPushButton("Start Docking Pipeline")
        self.btn_run.setStyleSheet("background-color: #2b8cbe; color: white; font-weight: bold; font-size: 14px; padding: 10px;")
        self.btn_run.clicked.connect(self.start_pipeline)
        
        self.btn_stop = QPushButton("Stop Screening")
        self.btn_stop.setEnabled(False) 
        self.btn_stop.clicked.connect(self.stop_pipeline)
        
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.btn_stop)
        
        self.progress_bar = QProgressBar()
        
        # Unified application stream console
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setPlaceholderText("Virtual screening processing steps and logs will stream here...")
        
        layout.addWidget(self.btn_run)
        layout.addLayout(controls_layout)
        layout.addWidget(QLabel("<b>Simulation Progress:</b>"))
        layout.addWidget(self.progress_bar)
        layout.addWidget(QLabel("<b>Pipeline Log Monitor:</b>"))
        layout.addWidget(self.log_output)
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "3. Run Simulation")

    def start_pipeline(self):
        rec = self.receptor_input.text()
        lig = self.ligand_input.text()
        if not rec or not lig:
            self.log_output.append("⚠️ Error: Select both receptor and ligand input fields.")
            return

        try:
            grid_config = {
                'center': [float(self.cx.text()), float(self.cy.text()), float(self.cz.text())],
                'size': [float(self.sx.text()), float(self.sy.text()), float(self.sz.text())],
                'exhaustiveness': int(self.txt_exhaustiveness.text()),
                'n_poses': int(self.txt_n_poses.text()),
                'cores_per_test_system': os.cpu_count() or 4
            }
        except ValueError:
            self.log_output.append("⚠️ Error: Check that all grid dimensions are numeric entries.")
            return

        self.log_output.clear()
        self.progress_bar.setValue(0)
        
        self.worker = DockingWorker(rec, lig, grid_config)
        self.worker.log_signal.connect(self.log_output.append)
        self.worker.progress_signal.connect(self.progress_bar.setValue)
        self.worker.finished_signal.connect(self.pipeline_finished)
        
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.worker.start()

    def stop_pipeline(self):
        if hasattr(self, 'worker') and self.worker.isRunning():
            self.worker.stop()
            self.btn_stop.setEnabled(False)
            self.log_output.append("\n🛑 Stop signal dispatched. Halting current calculations smoothly...")

    def pipeline_finished(self, final_msg):
        from PySide6.QtWidgets import QApplication
        self.log_output.append(f"\n{final_msg}")
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        QApplication.processEvents()
        self.load_results_table()
        self.tabs.setCurrentIndex(3)

    def process_pocket_calculation(self):
        if not self.chk_auto_center_pocket.isChecked(): return
        res_text = self.txt_residues.text().strip()
        if not res_text: return
        receptor_path = self.receptor_input.text()
        target_tokens = [t.strip().upper() for t in res_text.split(',') if t.strip()]
        x_pock, y_pock, z_pock = [], [], []
        try:
            with open(receptor_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM")):
                        res_name = line[17:20].strip().upper()
                        res_seq = line[22:26].strip()
                        if (f"{res_name}{res_seq}" in target_tokens) or (res_seq in target_tokens) or (res_name in target_tokens):
                            x_pock.append(float(line[30:38].strip()))
                            y_pock.append(float(line[38:46].strip()))
                            z_pock.append(float(line[46:54].strip()))
            if not x_pock: return
            self.cx.setText(f"{sum(x_pock)/len(x_pock):.3f}")
            self.cy.setText(f"{sum(y_pock)/len(y_pock):.3f}")
            self.cz.setText(f"{sum(z_pock)/len(z_pock):.3f}")
            self.update_volume_display()
        except Exception: pass

    def init_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel("<h2>Screening Performance Tables</h2>"))
        
        btn_refresh = QPushButton("🔄 Refresh Data")
        btn_refresh.clicked.connect(self.load_results_table)
        self.btn_export = QPushButton("📄 Export Data (.CSV)")
        self.btn_export.setStyleSheet("background-color: #2ca02c; color: white; font-weight: bold;")
        self.btn_export.clicked.connect(self.export_to_csv)
        
        top_layout.addWidget(btn_refresh)
        top_layout.addWidget(self.btn_export)
        layout.addLayout(top_layout)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Ligand Variant Name", "Binding Energy Profile (kcal/mol)"])
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.setSortingEnabled(True)
        layout.addWidget(self.table)

        self.btn_pymol = QPushButton("👁️ Deploy Pose Structure inside PyMOL Visualization Engine")
        self.btn_pymol.setStyleSheet("background-color: #00a86b; color: white; font-weight: bold; padding: 6px;")
        self.btn_pymol.clicked.connect(lambda: launch_pymol(self.table))
        layout.addWidget(self.btn_pymol)
        
        tab.setLayout(layout)
        self.tabs.addTab(tab, "4. Results Viewer")

    def load_results_table(self):
        output_dir = os.path.join(os.getcwd(), "docking_outputs")
        if not os.path.exists(output_dir): return
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0) 
        files = [f for f in os.listdir(output_dir) if f.lower().endswith('_results.pdbqt')]
        for file in files:
            ligand_name = file.replace('_results.pdbqt', '')
            affinity = "N/A"
            try:
                with open(os.path.join(output_dir, file), 'r') as f:
                    for line in f:
                        if "REMARK VINA RESULT:" in line:
                            affinity = float(line.split()[3])
                            break 
            except Exception: pass
            row_idx = self.table.rowCount()
            self.table.insertRow(row_idx)
            n_item = QTableWidgetItem(ligand_name); n_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            s_item = QTableWidgetItem()
            if isinstance(affinity, float): s_item.setData(Qt.DisplayRole, affinity)
            else: s_item.setText(affinity)
            s_item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.table.setItem(row_idx, 0, n_item)
            self.table.setItem(row_idx, 1, s_item)
        self.table.setSortingEnabled(True)
        self.table.sortItems(1, Qt.AscendingOrder)

    def export_to_csv(self):
        if self.table.rowCount() == 0: return
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Data Document", os.path.join(os.getcwd(), "screening_results.csv"), "CSV Files (*.csv)")
        if not file_path: return 
        try:
            with open(file_path, mode='w', newline='', encoding='utf-8') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["Ligand Name", "Best Binding Affinity (kcal/mol)"])
                for row in range(self.table.rowCount()):
                    writer.writerow([self.table.item(row, 0).text(), self.table.item(row, 1).text()])
            QMessageBox.information(self, "Export State", "File metrics written successfully.")
        except Exception as e: QMessageBox.critical(self, "System Crash Warning", str(e))

    def init_citation_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(QLabel("<h2>Citation Details</h2>"))
        citation_box = QTextEdit(); citation_box.setReadOnly(True)
        citation_text = """
        <h3>1. AutoDock Vina</h3>
        <p>Eberhardt, J., Santos-Martins, D., Tillack, A. F., & Forli, S. (2021). AutoDock Vina 1.2.0: New docking methods, expanded force field, and python bindings. <i>Journal of Chemical Information and Modeling</i>, 61(8), 3891-3898.</p>
        <h3>2. Open Babel</h3>
        <p>O'Boyle, N. M., Banck, M., James, C. A., Morley, C., Vandermeersch, T., & Hutchison, G. R. (2011). Open Babel: An open chemical toolbox. <i>Journal of Cheminformatics</i>, 3(1), 1-14.</p>
        """
        citation_box.setHtml(citation_text)
        layout.addWidget(citation_box)
        tab.setLayout(layout)
        self.tabs.addTab(tab, "Citation & References")