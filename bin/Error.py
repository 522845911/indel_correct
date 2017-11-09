#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True
class Error(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)
