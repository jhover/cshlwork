#!/usr/bin/env python
#
# consume CAFA prediction file. 
# 
#
#
#
#
  
import os
import sys
import logging

gitpath=os.path.expanduser("~/git/cshl-work")
sys.path.append(gitpath)

from expression.egad import *

print("cafaegad...")