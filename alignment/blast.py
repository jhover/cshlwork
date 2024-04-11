import logging
import os
import sys
import traceback

gitpath=os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

import numpy as np
import pandas as pd

from cshlwork.utils import run_command_shell, NonZeroReturnException, setup_logging
