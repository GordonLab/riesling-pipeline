"""
Clean, colorized logging for scripts.

Based on https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
with many adaptations by Nick Semenkovich <semenko@alum.mit.edu>

License: MIT
Author: Nick Semenkovich <semenko@alum.mit.edu>
"""

from __future__ import absolute_import, division, print_function, unicode_literals

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis'

import atexit
import datetime
import logging
import os
import platform
import sys
import time

if __name__ == '__main__':
    print("Do not run this as a standalone script.")
    exit()

# Very pretty error reporting, where available
try:
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Context', color_scheme='Linux')
except ImportError:
    pass

# Log our execution time, used by log_execution_time below.
STARTTIME = time.time()

## Define our core colors and resets.

RESET = "\x1b[0m"
BOLD = "\033[1m"
BLUE = "\x1b[34;01m"
CYAN = "\x1b[36;01m"
GREEN = "\x1b[32;01m"
RED = "\x1b[31;01m"
GRAY = "\x1b[37;01m"
YELLOW = "\x1b[33;01m"

# TODO: replace TMP with gettmp or whatever

def startLogger(verbose=False, noFileLog=False, initialLoggerName='main', outPath='/tmp'):
    """
    Set logging if called.
    TODO: Make this a class, extending the logging module.
    """
    datestamp = datetime.datetime.now().strftime("%Y%m%d-%H:%M:%S")
    results_path = outPath + '/' + datestamp + '.' + platform.node() + '.log'

    if not noFileLog:
        fileHandler = logging.FileHandler(results_path, 'w')
        fileFormatter = logging.Formatter('%(asctime)s: %(name)-25s: %(levelno)-3s: '
                                          '(%(filename)s:%(lineno)d, %(funcName)s) : %(message)s')
        fileHandler.setFormatter(fileFormatter)
        fileHandler.setLevel(logging.DEBUG)

    # Sneakily add colors directly as the log-level name.
    logging.addLevelName(logging.DEBUG,    CYAN   + 'DEBUG')
    logging.addLevelName(logging.INFO,     GREEN  + 'INFO')
    logging.addLevelName(logging.WARNING,  YELLOW + 'WARNING')
    logging.addLevelName(logging.ERROR,    RED    + 'ERROR')
    logging.addLevelName(logging.CRITICAL, RED    + 'CRITICAL')

    consoleHandler = logging.StreamHandler()
    consoleFormatter = logging.Formatter(BOLD + "%(name)-25s" + RESET + ": %(levelname)-17s" + RESET +
                                         ": %(message)-80s (" + BOLD + "%(filename)s" + RESET + ":%(lineno)d)")
    consoleHandler.setFormatter(consoleFormatter)

    consoleHandler.setLevel(logging.INFO)
    if verbose:
        consoleHandler.setLevel(logging.DEBUG)

    # Give out a logger!
    rootlog = logging.getLogger()
    rootlog.setLevel(logging.DEBUG)

    rootlog.addHandler(consoleHandler)
    if not noFileLog:
        rootlog.addHandler(fileHandler)
        rootlog.info('>> Logging to %s <<' % (results_path))
        rootlog.info('Running: %s [full command saved to log file]' % (os.path.basename(sys.argv[0])))
        rootlog.info('    Using python version: %s' % (sys.version.split('\n')[0]))
        rootlog.debug('%s' % (' '.join(sys.argv)))
        rootlog.info('Written by %s' % (__author__))
        rootlog.info('    Developed for the %s' % (__copyright__))

    return logging.getLogger(initialLoggerName)


def getLogger(name):
    """
    Returns a logger!
    """
    return logging.getLogger(name)

def log_execution_time():
    """
    Print the total elapsed execution time.

    :return:
    """
    logging.getLogger('root').info("Execution took: %0.2f secs." % (time.time() - STARTTIME))

# Register an exit handler to print execution time.
atexit.register(log_execution_time)
