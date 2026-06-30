import os
import subprocess
from PySide6.QtWidgets import QMessageBox



def launch_pymol(table_widget):
    """
    Extracts the selected ligand from a PySide6 QTableWidget,
    generates a custom .pml script, and launches it natively in PyMOL.
    """
    # 1. Extract selected text row coordinates
    selected_indexes = table_widget.selectedIndexes()
    if not selected_indexes:
        QMessageBox.warning(None, "No Selection", "Please click on a ligand row in the results table first.")
        return
        
    row = selected_indexes[0].row()
    # Assuming Column 0 contains the unique Ligand Name string identifier
    ligand_name = table_widget.item(row, 0).text() 
    
    output_dir = os.path.join(os.getcwd(), "docking_outputs")
    receptor_file = None
    
    # Verify the background path outputs exist
    if not os.path.exists(output_dir):
        QMessageBox.critical(None, "Directory Error", "The 'docking_outputs' directory does not exist.")
        return

    # Dynamic search for your prepared receptor file
    for f in os.listdir(output_dir):
        if f.endswith("_prepared.pdbqt"):
            receptor_file = os.path.join(output_dir, f)
            break
            


            
    ligand_file = os.path.join(output_dir, f"{ligand_name}_results.pdbqt")
    
    if not receptor_file or not os.path.exists(ligand_file):
        QMessageBox.critical(None, "Files Missing", f"Could not locate the structural coordinates for:\n{ligand_name}")
        return

    # 2. Write out a clean, automated PyMOL Script (.pml)
    script_path = os.path.join(output_dir, "view_interaction.pml")
    
    try:
        with open(script_path, "w", encoding="utf-8") as pml:
            pml.write("# Automated PyMOL Initialization Script\n")
            pml.write(f"load {receptor_file}, receptor\n")
            pml.write(f"load {ligand_file}, docked_ligand\n\n")
            
            # Format Protein Aesthetics (Clean grey cartoon ribbon)
            pml.write("hide everything, receptor\n")
            pml.write("show cartoon, receptor\n")
            pml.write("color gray80, receptor\n\n")
            
            # Format Ligand Aesthetics (High contrast magenta sticks)
            pml.write("hide everything, docked_ligand\n")
            pml.write("show sticks, docked_ligand\n")
            pml.write("color magenta, docked_ligand\n")
            pml.write("util.cnc docked_ligand\n\n") # Auto element colors (O=red, N=blue)
            
            # Highlight the local interaction pocket space (6.0 Angstrom view threshold)
            pml.write("select pocket, receptor within 6.0 of docked_ligand\n")
            pml.write("show sticks, pocket\n")
            pml.write("color gray40, pocket\n")
            pml.write("util.cnc pocket\n\n")
            
            # Highlight polar / hydrogen bonds
            pml.write("dist hbonds, docked_ligand, pocket, 3.5, mode=2\n")
            pml.write("hide labels, hbonds\n")
            pml.write("color yellow, hbonds\n\n")
            
            # Smooth camera centering directly onto the target zone
            pml.write("orient docked_ligand\n")
            pml.write("zoom docked_ligand, 8\n")
        
        # 3. Fire up the background native operating system subprocess call
        subprocess.Popen(["pymol", script_path])
        
    except Exception as e:
        QMessageBox.critical(None, "Visualizer Error", f"Failed to execute PyMOL sequence:\n{str(e)}")