#!/usr/bin/env python
from gi.repository import Gtk
import subprocess
import os

class dafuq:
	
	def __init__(self):
		path = os.path.dirname(os.path.realpath(__file__))
		filename = path+"/roda_gui.glade"
		builder = Gtk.Builder()
		builder.add_from_file(filename)
		builder.connect_signals(self)
		#self.text = builder.get_object("textbuffer1")
		self.window = builder.get_object("window1")
		self.window.connect("destroy", Gtk.main_quit)
		self.window.show()
		self.inputFile = builder.get_object("file_Input")
		self.outputFile = builder.get_object("file_Output")
		self.analysis = builder.get_object("combo_Analysis")
		self.analysis.set_active(0) # Glade fails to do this itself

	def analyse(self, widget):
		path = os.path.dirname(os.path.realpath(__file__))
		argInput=self.inputFile.get_filename()
		argOutput=self.outputFile.get_filename()
		argAnalysis=self.analysis.get_active_text()
		subprocess.Popen(args=["gnome-terminal", "-x", "bash", "-c", path+"/roda_mult.sh " + argInput + " " + argOutput + " " + argAnalysis + "; read -n1"]) 
		#subprocess.Popen(args=["gnome-terminal", "-x", "bash", "-c", "echo stuff; echo ; echo Press any key to close this window; read -n1"]) 

	def viewqueue(self, widget):
		command_view = "while true; do clear; qstat -s r -f -u '*'; echo ''; COUNT=`qstat -s p -u '*' | wc -l`; if [ ! $COUNT == 0 ]; then ((COUNT=$COUNT-2)); fi; echo 'Pending jobs:' $COUNT; sleep 3; done"
		subprocess.Popen(args=["gnome-terminal", "-x", "bash", "-c", command_view]) 

app = dafuq()
Gtk.main()
