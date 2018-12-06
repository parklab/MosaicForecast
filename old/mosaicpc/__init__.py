# -*- coding: utf-8 -*-
"""
MosaicPC
~~~~~~
:copyright: (c) 2016 Harvard Medical School
:author: Yanmei Dou, Minseok Kwon
:license: MIT
"""
import logging
__version__ = '0.7.11'
__format_version__ = 2
_logger = None

def get_logger():
    # Based on ipython traitlets
    global _logger

    if _logger is None:
        _logger = logging.getLogger('cooler')
        # Add a NullHandler to silence warnings about not being
        # initialized, per best practice for libraries.
        _logger.addHandler(logging.NullHandler())

    return _logger

