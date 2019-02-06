# -*- coding: utf-8 -*-
"""
Stack Overflow read eval print loop:
    
Created on Fri May 3 16:17:24 2010

@author: mathias
    http://stackoverflow.com/questions/2758159/how-to-embed-a-python-interpreter-in-a-pyqt-widget

when there is time to get this wired up:
    
    http://code.activestate.com/recipes/355319-using-codeinteractiveconsole-to-embed-a-python-she/

    https://docs.python.org/2/library/code.html

"""

import sys, os
import traceback
#
# did not work in 2014:
#from PyQt4 import QtCore
#from PyQt4 import QtGui
#import PyQt5
from PyQt5 import QtCore
from PyQt5 import QtGui
#if PyQt4:
#   QPlainTextEdit = QtGui.QPlainTextEdit
#elif PyQt5:
from PyQt5.QtWidgets import QPlainTextEdit
from PyQt5.QtWidgets import QApplication

#
# does not work in 2018:
#from PySide import QtCore, QtGui #does not work with iron python

#class Console(QtGui.QPlainTextEdit): #linux code (pyQt4)
class Console(QPlainTextEdit):
    def __init__(self, prompt='$> ', startup_message='', parent=None):
        #QtGui.QPlainTextEdit.__init__(self, parent) #linux code PyQt4
        QPlainTextEdit.__init__(self, parent)
        self.prompt = prompt
        self.history = []
        self.namespace = {}
        self.construct = []

        self.setGeometry(50, 60, 300, 100)
        self.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
        self.setUndoRedoEnabled(False)
        self.document().setDefaultFont(QtGui.QFont("monospace", 10, QtGui.QFont.Normal))
        self.showMessage(startup_message)
    
    
    def write_from_GUI(self, message):
        self.appendPlainText(message)
        #sys.stdout = stdoutProxy(message)
        self.newPrompt()
        #return 
        
    def keyPressEvent(self, e):
        """Only works for QtGui.QWidget"""
        if e.key() == QtCore.Qt.Key_Escape:
                self.close()
        return

    def updateNamespace(self, namespace):
        self.namespace.update(namespace)

    def showMessage(self, message):
        self.appendPlainText(message)
        self.newPrompt()

    def newPrompt(self):
        if self.construct:
            prompt = '.' * len(self.prompt)
        else:
            prompt = self.prompt
        self.appendPlainText(prompt)
        self.moveCursor(QtGui.QTextCursor.End)

    def getCommand(self):
        doc = self.document()
        curr_line = unicode(doc.findBlockByLineNumber(doc.lineCount() - 1).text())
        curr_line = curr_line.rstrip()
        curr_line = curr_line[len(self.prompt):]
        return curr_line

    def setCommand(self, command):
        if self.getCommand() == command:
            return
        self.moveCursor(QtGui.QTextCursor.End)
        self.moveCursor(QtGui.QTextCursor.StartOfLine, QtGui.QTextCursor.KeepAnchor)
        for i in range(len(self.prompt)):
            self.moveCursor(QtGui.QTextCursor.Right, QtGui.QTextCursor.KeepAnchor)
        self.textCursor().removeSelectedText()
        self.textCursor().insertText(command)
        self.moveCursor(QtGui.QTextCursor.End)

    def getConstruct(self, command):
        if self.construct:
            prev_command = self.construct[-1]
            self.construct.append(command)
            if not prev_command and not command:
                ret_val = '\n'.join(self.construct)
                self.construct = []
                return ret_val
            else:
                return ''
        else:
            if command and command[-1] == (':'):
                self.construct.append(command)
                return ''
            else:
                return command

    def getHistory(self):
        return self.history

    def setHisory(self, history):
        self.history = history

    def addToHistory(self, command):
        if command and (not self.history or self.history[-1] != command):
            self.history.append(command)
        self.history_index = len(self.history)

    def getPrevHistoryEntry(self):
        if self.history:
            self.history_index = max(0, self.history_index - 1)
            return self.history[self.history_index]
        return ''

    def getNextHistoryEntry(self):
        if self.history:
            hist_len = len(self.history)
            self.history_index = min(hist_len, self.history_index + 1)
            if self.history_index < hist_len:
                return self.history[self.history_index]
        return ''

    def getCursorPosition(self):
        return self.textCursor().columnNumber() - len(self.prompt)

    def setCursorPosition(self, position):
        self.moveCursor(QtGui.QTextCursor.StartOfLine)
        for i in range(len(self.prompt) + position):
            self.moveCursor(QtGui.QTextCursor.Right)

    def runCommand(self):
        command = self.getCommand()
        self.addToHistory(command)

        command = self.getConstruct(command)

        if command:
            tmp_stdout = sys.stdout

            class stdoutProxy():
                def __init__(self, write_func):
                    self.write_func = write_func
                    self.skip = False

                def write(self, text):
                    if not self.skip:
                        stripped_text = text.rstrip('\n')
                        self.write_func(stripped_text)
                        QtCore.QCoreApplication.processEvents()
                    self.skip = not self.skip

            sys.stdout = stdoutProxy(self.appendPlainText)
            try:
                try:
                    result = eval(command, self.namespace, self.namespace)
                    if result != None:
                        self.appendPlainText(repr(result))
                except SyntaxError:
                    exec command in self.namespace
            except SystemExit:
                self.close()
            except:
                traceback_lines = traceback.format_exc().split('\n')
                # Remove traceback mentioning this file, and a linebreak
                for i in (3,2,1,-1):
                    traceback_lines.pop(i)
                self.appendPlainText('\n'.join(traceback_lines))
            sys.stdout = tmp_stdout
        self.newPrompt()

    def keyPressEvent(self, event):
        if event.key() in (QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return):
            self.runCommand()
            return
        if event.key() == QtCore.Qt.Key_Home:
            self.setCursorPosition(0)
            return
        if event.key() == QtCore.Qt.Key_PageUp:
            return
        elif event.key() in (QtCore.Qt.Key_Left, QtCore.Qt.Key_Backspace):
            if self.getCursorPosition() == 0:
                return
        elif event.key() == QtCore.Qt.Key_Up:
            self.setCommand(self.getPrevHistoryEntry())
            return
        elif event.key() == QtCore.Qt.Key_Down:
            self.setCommand(self.getNextHistoryEntry())
            return
        elif event.key() == QtCore.Qt.Key_D and event.modifiers() == QtCore.Qt.ControlModifier:
            self.close()
        super(Console, self).keyPressEvent(event)

welcome_message = '''
   ---------------------------------------------------------------
     Welcome to a simple Python interpreter.
   ---------------------------------------------------------------
'''

if __name__ == '__main__':
    #app = QtGui.QApplication(sys.argv)
    app = QApplication(sys.argv)
    #console = Console(startup_message=welcome_message)
    console = Console()
    console.updateNamespace({'myVar1' : app, 'myVar2' : 1234})
    console.show();
    sys.exit(app.exec_())