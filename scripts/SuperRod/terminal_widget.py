import numpy as np
from pyqtgraph.console import ConsoleWidget
import pyqtgraph as pg

class TerminalWidget(ConsoleWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.localNamespace['pg'] = pg
        self.localNamespace['np'] = np

    def update_name_space(self, new_sub):
        self.localNamespace['model'] = new_sub

