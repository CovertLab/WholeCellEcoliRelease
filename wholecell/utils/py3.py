"""Python 2 -> 3 compatibility utilities beyond the six library."""

from __future__ import absolute_import, division, print_function

import time
from typing import Callable, Text, Union

__all__ = [
	"monotonic_seconds", "perf_counter_seconds", "process_time_seconds",
	"ANY_STRING", "String"]


#: This seconds clock for elapsed time relative measurements is unaffected by
#: system clock updates. The resolution might be nanoseconds. This might not be
#: seconds since the epoch, so `ctime()` might not format it correctly.
#: New in Python 3.3.
#: Fallback: `time.time()` reads the "Unix time" clock that can get set
#: backwards, doesn't count leap seconds, and has system-dependent precision
#: and meaning. The precision might not be better than 1 second.
monotonic_seconds = getattr(time, "monotonic", time.time)  # type: Callable[[], float]

#: This seconds clock for run time relative measurements has the highest
#: available resolution (probably nanoseconds) to measure short durations. It
#: includes time elapsed during sleep. This might not be seconds since the
#: epoch, so `ctime()` might not format it correctly.
#: New in Python 3.3.
#: Fallback: `time.time()`.
perf_counter_seconds = getattr(time, "perf_counter", time.time)  # type: Callable[[], float]

#: This seconds clock for CPU time relative measurements includes system + user
#: CPU time for the current process and excludes sleep time. This might not be
#: seconds since the epoch, so `ctime()` might not format it correctly.
#: New in Python 3.3.
#: Fallback: `time.time()`.
process_time_seconds = getattr(time, "process_time", time.time)  # type: Callable[[], float]

# time.clock()  -- This clock's precision and meaning are system-dependent, its
# precision might not be better than 1 second, and it can get set backwards.
# It was deprecated in Python 3.3 and removed in 3.8.


#: Any type of string, for use with `isinstance()`.
#: After dropping Python 2, this can become `(bytes, str)`.
ANY_STRING = (bytes, str, Text)

#: A type alias for Python 2 str or unicode; Python 3 str (not bytes), for use
#: in type hints. [`isinstance()` doesn't accept Unions so use a tuple there
#: like `ANY_STRING` or `six.string_types`.]
#: After dropping Python 2, this can become `str`.
String = Union[str, Text]
