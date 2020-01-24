from PyQt5 import QtWidgets
import syntax_pars

app = QtWidgets.QApplication([])
editor = QtWidgets.QPlainTextEdit()
editor.setStyleSheet("""QPlainTextEdit{
	font-family:'Consolas'; 
	color: #ccc; 
	background-color: #2b2b2b;}""")
highlight = syntax_pars.PythonHighlighter(editor.document())
editor.show()

# Load syntax.py into the editor for demo purposes
# infile = open('C:\\apps\\DaFy_P23\scripts\\SuperRod\\syntax_pars.py', 'r')
# editor.setPlainText(infile.read())
editor.setPlainText('')

app.exec_()