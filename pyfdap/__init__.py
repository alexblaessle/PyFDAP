#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFDAP: automated analysis of
Fluorescence Decay After Photoconversion
(FDAP) experiments
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time
import csv
import pickle
import sys
import platform

from . import pyfdap_misc_module
from . import pyfdap_img_module
from . import pyfdap_fit_module
from . import pyfdap_stats_module
from . import embryo
from . import molecule
from . import useful_fcts

from . import pyfdap_app
from . import pyfdap_plot_dialogs

__version__ = '1.1.1'
__author__ = u"Alexander Blaessle"
__license__ = "GNU GPL v3"
