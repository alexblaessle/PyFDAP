#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller, and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFDAP.

#PyFDAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.

#=====================================================================================================================================
#Module Description
#=====================================================================================================================================

#Module containing classes to emulate Python terminal in PyQT GUI for PyFRAP.
#Modified from http://stackoverflow.com/questions/2758159/how-to-embed-a-python-interpreter-in-a-pyqt-widget/8219536#8219536

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

import os
import re
import sys
import code
import warnings

from PyQt4.QtGui import *
from PyQt4.QtCore import *

#=====================================================================================================================================
#Module classes
#=====================================================================================================================================

class PyInterp(QTextEdit):

    class InteractiveInterpreter(code.InteractiveInterpreter):

        def __init__(self, locals):
            code.InteractiveInterpreter.__init__(self, locals)
            
            
        def runIt(self, command,ignWarnings=True):

            if ignWarnings:
	        with warnings.catch_warnings():
		   
		    warnings.filterwarnings("ignore",category=DeprecationWarning)
		    warnings.filterwarnings("ignore",category=FutureWarning)	
		
		    code.InteractiveInterpreter.runsource(self, command)
            else:
		    code.InteractiveInterpreter.runsource(self, command)
                    

    def __init__(self,  parent,ignWarnings=True,redirect=True):
        super(PyInterp,  self).__init__(parent)
	
	if redirect:
	     sys.stdout              = self
	     sys.stderr              = self
        self.refreshMarker      = False # to change back to >>> from ...
        self.multiLine          = False # code spans more than one line
        self.command            = ''    # command to be ran
        self.ignWarnings = ignWarnings
        
        self.history            = []    # list of commands entered
        self.historyIndex       = -1
        self.interpreterLocals  = {}
	
        # setting the color for bg and text
        palette = QPalette()
        palette.setColor(QPalette.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.Text, QColor(0, 0, 0))
        self.setPalette(palette)
        self.setFont(QFont('Courier', 12))

        # initilize interpreter with self locals
        self.initInterpreter(locals())
	
	self.interpreter.runIt('import pyfdap_img_module as img',ignWarnings=self.ignWarnings)
	self.interpreter.runIt('import pyfdap_fit_module as fit',ignWarnings=self.ignWarnings)
	self.interpreter.runIt('import pyfdap_misc_module as misc',ignWarnings=self.ignWarnings)
	self.interpreter.runIt('from numpy import *',ignWarnings=self.ignWarnings)
	self.interpreter.runIt('import useful_fcts as uf',ignWarnings=self.ignWarnings)
	
		
	self.printBanner()    
	self.marker()                   # make the >>> or ... marker     
	
    def printBanner(self):
        self.write(sys.version)
        self.write(' on ' + sys.platform + '\n')
        self.write('PyQt4 ' + PYQT_VERSION_STR + '\n')
          
        self.write("---------------------------------------"+ '\n')
        msg = 'Type !hist for a history view and !hist(n) history index recall'
        self.write(msg + '\n')
        msg2 = 'To access local variables, type pyfdp.variable. For example, to access the currently selected molecule, type pyfdp.curr_mol. For more details, see the documentation.'
        self.write(msg2 + '\n')
        msg3 = 'To access functions of the pyfdp modules, type modulename.function(). Available modulenames are "fit", "img" and "misc".'
        self.write(msg3 + '\n')
        self.write("---------------------------------------" + '\n')

    def marker(self):
        if self.multiLine:
            self.insertPlainText('... ')
        else:
            self.insertPlainText('>>> ')

    def initInterpreter(self, interpreterLocals=None):
     
        if interpreterLocals:
            # when we pass in locals, we don't want it to be named "self"
            # so we rename it with the name of the class that did the passing
            # and reinsert the locals back into the interpreter dictionary
            
            selfName = interpreterLocals['self'].__class__.__name__
            interpreterLocalVars = interpreterLocals.pop('self')
            self.interpreterLocals[selfName] = interpreterLocalVars
                      
        else:
            self.interpreterLocals = interpreterLocals
        self.interpreter = self.InteractiveInterpreter(self.interpreterLocals)

    def updateInterpreterLocals(self, newLocals):
        className = newLocals.__class__.__name__
        self.interpreterLocals[className] = newLocals

    def write(self, line):
        self.insertPlainText(line)
        self.ensureCursorVisible()

    def clearCurrentBlock(self):
        # block being current row
        length = len(self.document().lastBlock().text()[4:])
        if length == 0:
            return None
        else:
            # should have a better way of doing this but I can't find it
            [self.textCursor().deletePreviousChar() for x in xrange(length)]
        return True

    def recallHistory(self):
        # used when using the arrow keys to scroll through history
        self.clearCurrentBlock()
        if self.historyIndex <> -1:
            self.insertPlainText(self.history[self.historyIndex])
        return True

    def customCommands(self, command):

        if command == '!hist': # display history
            self.append('') # move down one line
            # vars that are in the command are prefixed with ____CC and deleted
            # once the command is done so they don't show up in dir()
            backup = self.interpreterLocals.copy()
            history = self.history[:]
            history.reverse()
            for i, x in enumerate(history):
                iSize = len(str(i))
                delta = len(str(len(history))) - iSize
                line = line  = ' ' * delta + '%i: %s' % (i, x) + '\n'
                self.write(line)
            self.updateInterpreterLocals(backup)
            self.marker()
            return True

        if re.match('!hist\(\d+\)', command): # recall command from history
            backup = self.interpreterLocals.copy()
            history = self.history[:]
            history.reverse()
            index = int(command[6:-1])
            self.clearCurrentBlock()
            command = history[index]
            if command[-1] == ':':
                self.multiLine = True
            self.write(command)
            self.updateInterpreterLocals(backup)
            return True

        return False

    def keyPressEvent(self, event):

        if event.key() == Qt.Key_Escape:
            # proper exit
            self.interpreter.runIt('exit()',ignWarnings=self.ignWarnings)

        if event.key() == Qt.Key_Down:
            if self.historyIndex == len(self.history):
                self.historyIndex -= 1
            try:
                if self.historyIndex > -1:
                    self.historyIndex -= 1
                    self.recallHistory()
                else:
                    self.clearCurrentBlock()
            except:
                pass
            return None

        if event.key() == Qt.Key_Up:
            try:
                if len(self.history) - 1 > self.historyIndex:
                    self.historyIndex += 1
                    self.recallHistory()
                else:
                    self.historyIndex = len(self.history)
            except:
                pass
            return None

        if event.key() == Qt.Key_Home:
            # set cursor to position 4 in current block. 4 because that's where
            # the marker stops
            blockLength = len(self.document().lastBlock().text()[4:])
            lineLength  = len(self.document().toPlainText())
            position = lineLength - blockLength
            textCursor  = self.textCursor()
            textCursor.setPosition(position)
            self.setTextCursor(textCursor)
            return None

        if event.key() in [Qt.Key_Left, Qt.Key_Backspace]:
            # don't allow deletion of marker
            if self.textCursor().positionInBlock() == 4:
                return None

        if event.key() in [Qt.Key_Return, Qt.Key_Enter]:
            # set cursor to end of line to avoid line splitting
            textCursor = self.textCursor()
            position   = len(self.document().toPlainText())
            textCursor.setPosition(position)
            self.setTextCursor(textCursor)

            line = str(self.document().lastBlock().text())[4:] # remove marker
            line.rstrip()
            self.historyIndex = -1

            if self.customCommands(line):
                return None
            else:
                try:
                    line[-1]
                    self.haveLine = True
                    if line[-1] == ':':
                        self.multiLine = True
                    self.history.insert(0, line)
                except:
                    self.haveLine = False

                if self.haveLine and self.multiLine: # multi line command
                    self.command += line + '\n' # + command and line
                    self.append('') # move down one line
                    self.marker() # handle marker style
                    return None

                if self.haveLine and not self.multiLine: # one line command
                    self.command = line # line is the command
                    self.append('') # move down one line
                    self.interpreter.runIt(self.command,ignWarnings=self.ignWarnings)
                    self.command = '' # clear command
                    self.marker() # handle marker style
                    return None

                if self.multiLine and not self.haveLine: #  multi line done
                    self.append('') # move down one line
                    self.interpreter.runIt(self.command,ignWarnings=self.ignWarnings)
                    self.command = '' # clear command
                    self.multiLine = False # back to single line
                    self.marker() # handle marker style
                    return None

                if not self.haveLine and not self.multiLine: # just enter
                    self.append('')
                    self.marker()
                    return None
                return None

        # allow all other key events
        super(PyInterp, self).keyPressEvent(event)
