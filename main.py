import sys
from PySide6.QtWidgets import QApplication
from gui import DockingApp

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DockingApp()
    window.show()
    sys.exit(app.exec())


    