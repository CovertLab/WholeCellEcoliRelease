# Style Guide

_These guidelines are condensed from the PEP8 and Google style guides here, along with our decisions. -- Jerry_

* See [PEP8](https://www.python.org/dev/peps/pep-0008/)
* See [Google's Python style guide](https://google.github.io/styleguide/pyguide.html)
* See our [March 12, 2018 meeting slides](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_55) with notes on the tradeoffs between alternatives.

## Note: Python 2 and 3 compatibility
Since Python 2 end-of-life is Jan. 1, 2020, the project is transitioning to Python 3. **Code additions and changes should remain compatible with Python 2 and also be compatible with Python 3 to an increasing degree.**

* Until we finish the transition, use syntax and features that are compatible with both Python 2 + 3. We'll use the conversion tools to do the bulk of the conversion work, but do make edits compatible. Easy cases:
  *  `raise Exception, args` → `raise Exception(args)`
  *  `except TypeError, ZeroDivisionError:` → `except Exception as variable:`
  *  `` `i` `` → `repr(i)`
  *  `a <> b` → `a != b`
  *  `a_dict.has_key(k)` → `k in a_dict`
  *  `apply(f, args)` → `f(*args)`

* Put `from __future__ import absolute_import, division, print_function` in each source file after making and testing any needed adaptations.
* See [Python 2.7 → 3.7 slides](https://docs.google.com/presentation/d/1xwVHjpQsRTGLdaIpPDNjKPQFtp7mBNfn7oTWK4Kct30/edit?usp=sharing) and [TODO] wiki pages.


# Style Guides

Style guides make recommendations among programming alternatives like [imports](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_107), [docstrings](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_123), [common patterns](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_117), [names](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_151), and [layout](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_131). The purpose is to reduce the likelihood of some bugs, increase code readability and familiarity, and make it easier for programmers to collaborate and merge code. But don't overdo consistency.

For each guideline, we could decide to:
1. Set an expectation. Point it out in code reviews.
2. Soft target. Don't sweat it in code reviews.
3. Don't adopt it.

and we can plan to:
1. Adopt it for new code.
2. Plan to change existing code gradually or rapidly.


## PEP8 Style Guidelines -- with a few adjustments

* **Indentation:** Stick with TABs in this project.
  * Set your editor's TAB stops to 4 spaces.
  * Python assumes TAB stops are at 8 spaces unless told otherwise.
  * Python 3 disallows mixing the use of tab and space indentation.
  * Set a nightly build script or a check-in script to flag SPACE indentation. It could use `python -t sourcefile.py`, or `tabnanny` which is less lenient but still allows mixtures that Python allows, or just search for any SPACE in indentation (although it's normal to use TABs followed by some SPACEs, esp. for odd indentation to line up with a `(` or for half-TAB indentation when tab stops are 8 spaces):

        find . -name "*.py" -o -name "*.pyx" | xargs grep '^\s* '

    although this will include indentation in multi-line comments and strings and indentation that aligns a continuation line with the previous line's `(` or `[`.

  * (PEP8 recommends 4 spaces per indentation level but that's primarily for _shared code_.)

* Use 7-bit ASCII text in Python 2. It will be interpreted as valid UTF-8 in Python 3. Otherwise the file needs an encoding declaration.

* The **line length** soft target is 79 columns; harder target at 99 columns; no hard limit. The same limit for comments and docstrings.
  * A standard line length aids editor window layout and diff displays, but bio simulations might have many long names. It's annoying to stay within a hard limit but very useful to have a shared target.
  * (PEP8 recommends 79 columns, but 72 for comments and docstrings.)
  * A shell script to check for very long lines in source files:

        find . -name '*.py' -exec awk '{ if (length($0) > max) max = length($0) } END { if (max > 199) print max, FILENAME }' {} \;

* Don't use implicit relative imports (e.g. `import sibling` where `sibling` is in the same directory) because it can import the wrong file (e.g. `import random`), it can import the same module twice (really?), and it doesn't work in Python 3.

  Instead use absolute imports or explicit relative imports:

      from __future__ import absolute_import  # prevents implicit relative imports
      from . import sibling
      from path.to.mypkg import sibling
      from .sibling import example

* Put **imports at the top of the file.** Certainly don't do `import` in repeatedly-executed code.

  Occasionally there are good reasons to break this rule, like `import pdb`.

* Import separate modules on separate lines.

* **Avoid wildcard imports** (`from <module> import *`).
  * Never `import *` within a class or a function. That generates slow code and it won't compile in Python 3.

* Use `if x is None:` or `if x is not None:` rather than `==` or `!=`, and likewise for other singletons like enum values (see pip enum34). It states a simpler aim. It's faster, esp. if it avoids calling a custom `__eq__()` method, and it might avoid exceptions or incorrect results in `__eq__(None)`.

* Prefer to use Python's **implied line continuation** inside parentheses, brackets and braces over a backslash for line continuation.

* Write **docstrings** for all public modules, functions, classes, and methods.

  A function docstring should start with an imperative sentence like "Plot all the things.".

  Add docstrings to existing classes and functions while working on the code. A short introduction can be a big help vs. nothing.

  Comments and docstrings that contradict the code are worse than no comments.

* Reevaluate what is `public` vs. `_private`.


* **Line continuation**

  preferred:

      def long_function_name(
              var_one, var_two, var_three,
              var_four):
          print(var_one)

  preferred:

      def long_function_name(
              var_one, var_two, var_three,
              var_four
              ):
          print(var_one)

  PEP8 alternative:

      foo = long_function_name(var_one, var_two,
                               var_three, var_four)

  PyCharm is configurable but it implements the first indentation style by default, and using its Refactor command to rename "long_function_name" will automatically adjust the indentation of the continuation lines.

* Import order: standard library imports, blank line, third party imports, blank line, local imports. It's nice to alphabetize each group by module name.

* Module content order:

      """Module docstring."""
      from __future__ import ...
      __all__ = ['a', 'b', 'c']  # and other "dunder" settings like __version__ and __author__
      imports
      code

* `"""Docstrings in triple quotes."""` whether it's one line or many; whether `"""` or `'''`. The `"""` that ends a multiline docstring should be on a line by itself.

* Two blank lines between top-level definitions. One blank line between method definitions. An additional blank line is OK to separate sections.

* Prefer to put a line break before a binary operator, but after is also OK.

* Put at least one space before an inline comment, then `#␣` (that's one space after the `#`). (PEP8 says at least two spaces before the `#`, but use your judgement.)

* Spacing like this (see [PEP8](https://www.python.org/dev/peps/pep-0008/) for more info):

      # Put no spaces immediately within `()`, `[]`, or `{}`.
      spam(ham[1], {eggs: 2, salsa: 10})
      
      # Put a space between `,` `;` or `:` and any following item.
      demo = (0,) + (2, 3)
      if x == 4:
          print x, y; x, y = y, x
      
      # Put no space in a simple slice expression, but parentheses to clarify complicated slice precedence
      # or construct a slice object or put subexpressions into variables.
      ham[1:9], ham[1:9:3], ham[:9:3], ham[1::3], ham[1:9:]
      ham[(lower+offset):(upper+offset)]
      
      # Put no space in function application or object indexing.
      spam(1) < spam(2)
      dct['key'] += lst[index]
      
      # Don't line up the `=` on multiple lines of assignment statements.
      x = 1
      long_variable = (3, 10)
      
      # Spaces around keyword `=` are OK, unlike in PEP8, which recommends them only
      # when there's a Python 3 parameter annotation.
      c = magic(real=1.0, imag=10.5)
      c = magic(real = 1.0, imag = 10.5)
      def munge(input: AnyStr, sep: AnyStr = None, limit=1000): ...
      
      # Use spaces or parentheses to help convey precedence.
      # Put zero or one space on both sides of a binary operator (except indentation).
      hypot2 = x*x + y*y

  Avoid trailing whitespace -- a backslash followed by a space and a newline does not count as a line continuation marker.

* Avoid compound statements on one line.

      if foo == 'blah': do_something()

* Comments are usually complete sentences. The first word should be capitalized unless it's an identifier that begins with a lower case letter.

* Stylize **names** like this:

   ```python
   ClassName
   ExceptionName  # usually ends with "Error"
   
   GLOBAL_CONSTANT_NAME
   
   function_name, method_name
   decorator_name
   
   local_var_name, global_var_name, instance_var_name, function_parameter_name
   camelCase  # OK to match the existing style
   
   __mangled_class_attribute_name
   _internal_name
   
   module_name
   package  # underscores are discouraged
   ```

    * Public names (like a class used as a decorator) follow conventions for usage rather than implementation.
    * Use a trailing `_` to avoid conflicting with a Python keyword like `yield_`, `complex_`, or `max_`.
    * Don't invent `__double_leading_and_trailing_underscore__` special names.
    * Always use `self` for the first argument of an instance method and `cls` for the first argument of a class method.
    * Don't use `l`, `O`, or `I` for single character variable names.
    * Don't make exceptions for scientific conventions like `Kcat` and math conventions like matrix `M`, and any name is better than a single letter.
    * Avoid using properties for expensive operations. The attribute notation suggests it's cheap.
    * Use the verb to distinguish methods like `get_value()` from `compute_value()`.

* Documented interfaces are considered public unless the documentation says they're provisional or internal. Undocumented interfaces are assumed to be internal.

* The `__all__` attribute is useful for introspection.

Programming tips:

* `if x is not None` is more readable than `if not x is None`.
* When implementing ordering operations with rich comparisons, it's best to implement all six operations or use the `functools.total_ordering()` decorator to fill them out.
* Use `def f(x): return 2*x` instead of `f = lambda x: 2*x` for more helpful stack traces.
* Derive exceptions from Exception rather than BaseException unless catching it is almost always the wrong thing to do.
* When designing and raising exceptions aim to answer the question "What went wrong?" rather than only indicating "A problem occurred."
* In Python 2, use `raise ValueError('message')` instead of `raise ValueError, 'message'` (which is not legal in Python 3).
* Use the bare `except:` clause only when printing/logging the traceback.
* Use the form `except Exception as exc:` to bind the exception name.
* Limit a `try` clause to a narrow range of code so it only doesn't bury totally unexpected exceptions.
* Use a `with` statement or try/finally to ensure cleanup gets done. For a file-like object that that doesn't support the `with` statement, use `with contextlib.closing(urllib.urlopen("https://www.python.org/")):`.
* In a function, make either all or none of the `return` statements return an explicit value.
  * Furthermore, have a consistent return type. Make a class instance, tuple, `namedtuple`, or dictionary to handle a union of different cases.
  * Any kind of failure should raise an explicit exception.

* Use string methods instead of the string module. They're faster and have the same API as Unicode strings.
* String `.startswith()` and `.endswith()` are less error prone than string slicing.
* Use e.g. `isinstance(obj, int)` instead of `type(obj) is type(1)` to check an object's type. Use `isinstance(obj, basestring)` to accept both str and unicode.
  * Better yet, avoid checking types except to catch common errors. It's cleaner to call different function for distinct input patterns or use O-O dispatch.

* Use `' '.join()` rather than looping over `a_string += stuff` to combine strings since `join()` takes linear time rather than O(_n_^2) time.



## From Google's Python Style Guide

See https://google.github.io/styleguide/pyguide.html

* Use pylint. [And/or PyCharm inspections.]
* Import packages and modules only, not names from modules.
* Use full pathnames to import a module; no relative imports to help prevent importing a package twice.
* Avoid global variables.
* Use the operator module e.g. `operator.mul` over `lambda x, y: x * y`. There's also `operator.itemgetter(*items)`, `operator.attrgetter(*attrs)`, and `operator.methodcaller(name[, args...])`.
* Don't use mutable objects as default values in a function or method definition.
* Use Python falsy tests for empty sequences and 0, e.g. `if sequence:` rather than `if len(sequence):` or `if len(sequence) > 0:`, but not for testing if a value is (not) `None`.
  * But don't write `if value:` to test for a non-empty string. That can be confusing.
* Avoid (or minimize) features such as metaclasses, access to bytecode, on-the-fly compilation, dynamic inheritance, object reparenting, import hacks, reflection, modification of system internals, etc. [at least without compelling reasons].
* Use extra parentheses instead of backslash line continuation.
* Don't use parentheses in return statements or conditional statements except for implied line continuation or tuples.
  * Prefer explicit tuple parentheses, definitely for 1-element tuples. `(x,) = func()` not `x, = func()`. `for (i, x) in enumerate(...):`.
  * PyCharm recommends `return x, y` over `return (x, y)`.
  * There's a case to make for always using parentheses for tuple construction and tuple unpacking, but we haven't set a guideline for that.

* Write a docstring as a summary sentence, a blank line, then the rest.
* A function must have a docstring, unless it's: not externally visible, short, and obvious. Explain what you need to know to call it.
* Classes should have docstrings.
* Don't put code at the top level that you don't want to run when the file is imported, unit tested, or pydoc'd.



## Other Recommendations

* Prefer to write floating point literals like `1.0` or `0.1` for clarity rather than `1.` or `.1`, but NumPy uses `1.` and `0.1` so we're not ruling that out.

* A function should return consistent types rather than different types in different cases. Construct a dict or object or tuple to handle a union of cases. This makes it easier to understand and use.

* Prefer `format_string.format(...)` over `format_string % ...` because [printf-style % formatting has](https://docs.python.org/3/library/stdtypes.html?highlight=sprintf#printf-style-string-formatting) "a variety of quirks that lead to a number of common errors (such as failing to display tuples and dictionaries correctly)". `.format()` is a bit more readable and it offers some nice additional features. (See: More info on [the 4 ways to format strings in Python](https://dbader.org/blog/python-string-formatting).)

# Performance Tips

[Summarized from sources like [PythonSpeed](https://wiki.python.org/moin/PythonSpeed).]

* Testing membership in a set or a dict is very fast, `O(n)`, unlike a list, tuple, or array.
* Sorting a list using a _sort key_ is faster than using a _comparison function._
* Mapping a function over a list, or using a list comprehension or generator comprehension, should be faster than a `for` loop since it pushes the loop work into compiled C code.
* Local variables are faster to access than global variables, builtins, and attribute lookups.
* Iterators are generally more memory-friendly and scalable than list operations, e.g. `xrange()`, `itertools.imap()`, `dict.iteritems` vs. `range()`, `map()`, `dict.items()`.
* Core building blocks are coded in optimized C, including the builtin datatypes (lists, tuples, sets, and dictionaries) and extension modules like `array`, `itertools`, and `collections.deque`.
* Builtin functions run faster than hand-built equivalents, e.g. `map(operator.add, v1, v2)` is faster than `map(lambda x, y: x+y, v1, v2)`.
* For queue applications using `pop(0)` or `insert(0,v)`, `collections.deque()` offers superior `O(1)` performance over a list because it avoids the `O(n)` step of rebuilding a list for each insertion or deletion.
* Chained comparisons like `x < y < z` are faster and hopefully more readable than `x < y and y < z`.
* Threading can improve the response time in applications that would otherwise waste time waiting for I/O.
