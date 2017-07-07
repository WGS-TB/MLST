# --------------------------------------------------------------------------
# File: _ostream.py 
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2015. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------
"""
"""

import weakref

from ._procedural import check_status
from ..exceptions import CplexError
from .. import six

class OutputStream(object):
    """Class to parse and write strings to a file object.

    Methods:
    __init__(self, outputfile, fn = None)
    __del__(self)
    write(self)
    flush(self)
    """

    def __init__(self, outputfile, env, fn=None, initerrorstr=False):
        """OutputStream constructor.

        outputfile must provide methods write(self, str) and
        flush(self).

        If fn is specified, it must be a fuction with signature
        fn(str) -> str.
        """
        self._env = weakref.proxy(env)
        self._fn = fn
        self._is_valid = False
        self._was_opened = False
        self._disposed = False
        # We only create this attribute for the error channel.
        if initerrorstr:
            self._error_string = None
        if isinstance(outputfile, six.string_types):
            self._file = open(outputfile, "w")
            self._was_opened = True
        else:
            self._file = outputfile
        if self._file is not None:
            try:
                tst = callable(self._file.write)
            except AttributeError:
                tst = False
            if not tst:
                raise CplexError("Output object must have write method")
            try:
                tst = callable(self._file.flush)
            except AttributeError:
                tst = False
            if not tst:
                raise CplexError("Output object must have flush method")
        self._is_valid = True

    def _end(self):
        """Flush and free any open file.

        If the user passes in a filename string, we open it.  In that case,
        we need to clean it up.
        """
        if self._disposed:
            return
        # File-like objects should implement this attribute.  If we come
        # across one that doesn't, don't assume anything.
        try:
            isclosed = self._file.closed
        except AttributeError:
            isclosed = False
        # If something bad happened in the constructor, then don't flush.
        if self._is_valid and not isclosed:
            self.flush()
            # If we opened the file, then we need to close it.
            if self._was_opened:
                self._file.close()
        self._disposed = True

    def __del__(self):
        """OutputStream destructor."""
        self._end()

    def _write_wrap(self, str_):
        """Only used by callbacks (see SWIG_callback.c)."""
        try:
            self._terminate = 0
            self.write(str_)
            self.flush()
            try:
                msg = self._error_string
            except AttributeError:
                msg = None
            if msg is not None:
                if not msg.startswith("CPLEX Error  1006"):
                    self._error_string = None
                    raise CplexError("ERROR", msg)
        except Exception as exc:
            self._env._callback_exception = exc
            check_status._pyenv = self._env
            self._terminate = 1

    def write(self, str_):
        """Parses and writes a string.

        If self._fn is not None, self._fn(str_) is passed to
        self._file.write.  Otherwise, str_ is passed to self._file.write
        """
        if self._file is None:
            return
        if str_ is None:
            str_ = ""
        if self._fn is None:
            self._file.write(str_)
        else:
            self._file.write(self._fn(str_))

    def flush(self):
        """Flushes the buffer."""
        if self._file is not None:
            self._file.flush()
