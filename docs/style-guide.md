# Style Guide

_These guidelines are derived from the PEP8 and Google style guides._

* See [PEP8](https://www.python.org/dev/peps/pep-0008/)
* See [Google's Python style guide](https://google.github.io/styleguide/pyguide.html)

_NOTE: This release snapshot was partway through the project's transition to Python 3.
This code runs in Python 2._


Style guides make recommendations among programming alternatives like [imports](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_107), [docstrings](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_123), [common patterns](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_117), [names](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_151), and [layout](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_131). The purpose is to reduce the likelihood of some bugs, increase code readability and familiarity, and make it easier for programmers to collaborate and merge code. But don't overdo consistency.


## PEP8 Style Guidelines -- with a few adjustments

* **Indentation:** Stick with TABs in this project.
  * Set your editor's TAB stops to 4 spaces.
  * Python assumes TAB stops are at 8 spaces unless told otherwise.
  * Python 3 disallows mixing the use of tab and space indentation.
  * (PEP8 recommends 4 spaces per indentation level but that's primarily for _shared code_.)

* Use 7-bit ASCII text in Python 2. It will be interpreted as valid UTF-8 in Python 3. Otherwise a file needs an encoding declaration.

* The **line length** soft target is 79 columns; harder target at 99 columns; no hard limit. The same limit for comments and docstrings.
  * A standard line length aids editor window layout and diff displays, but bio simulations might have many long names. It's annoying to stay within a hard limit but very useful to have a shared target.
  * (PEP8 recommends 79 columns, but 72 for comments and docstrings.)

* Don't use implicit relative imports (e.g. `import sibling` where `sibling` is in the same directory) because it can import the wrong module (e.g. `import random`) and it doesn't work in Python 3.

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

* Use `if x is None:` or `if x is not None:` rather than `==` or `!=` for `None` and other singletons like enum values (see pip enum34). It states a simpler aim. It's faster, esp. if it avoids calling a custom `__eq__()` method, and it might avoid exceptions or incorrect results in `__eq__(None)`.

* Prefer to use Python's **implied line continuation** inside parentheses, brackets, and braces over a backslash for line continuation.

* Write **docstrings** for all public modules, functions, classes, and methods.

  A function docstring should start with an imperative sentence like "Plot all the things.".

  Add docstrings to existing classes and functions while working on the code. Even a short introduction can be a big help.

  Comments and docstrings that contradict the code are worse than no comments.

* Reevaluate what to make `public` vs. `_private`.


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

* Put at least one space before an inline comment, preferably two spaces, then `#␣` (that's one space after the `#`).

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
      
      # Prefer no spaces around keyword_parameter=default. PEP8 recommends such spaces only
      # when there's a Python 3 parameter annotation.
      c = magic(real=1.0, imag=10.5)
      c = magic(real = 1.0, imag = 10.5)
      def munge(input: AnyStr, sep: AnyStr = None, limit=1000): ...
      
      # It's OK to use spaces or extra parentheses to help convey precedence.
      hypot2 = x*x + y*y

  Avoid trailing whitespace -- a backslash followed by a space and a newline does not count as a line continuation marker.

* Avoid compound statements on one line.

      if foo == 'blah': do_something()

* Comments are usually complete sentences. The first word should be capitalized unless it's an identifier that begins with a lower case letter.

  Per PEP 257, write docstrings in the command form: `"Do this."` `"Return that."`

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
* Use Python 3 compatible syntax: `raise ValueError('message')` and `except Exception as exc:`.
* Use the bare `except:` clause only when printing/logging the traceback.
* Catch a narrow set of exceptions rather than unintentionally catch a
  `KeyboardInterrupt`, `SystemExit`, `NameError`, etc.
* Limit a `try` clause to a narrow range of code so it only doesn't bury totally unexpected exceptions.
* Use a `with` statement or `try`/`finally` to ensure cleanup gets done. For a file-like object that doesn't support the `with` statement, use `contextlib.closing`,
  e.g. `with contextlib.closing(urllib.urlopen("https://www.python.org/")):`.
* In a function, make either all or none of the `return` statements return an explicit value.
  * A function should return consistent types rather than different types in different cases. Construct a class instance, tuple, `namedtuple`, or dictionary as needed.
  * A failure case should raise an explicit exception.
* Use string methods instead of the string module. They're faster and have the same API as Unicode strings.
* String `.startswith()` and `.endswith()` are less error-prone than string slicing.
* Use e.g. `isinstance(obj, int)` instead of `type(obj) is type(1)` to check an object's type. Use `isinstance(obj, basestring)` to accept both str and unicode.
  * Better yet, avoid checking types except to catch common errors. It's cleaner to call different function for distinct input patterns or use O-O dispatch.
* Use `' '.join()` rather than looping over `a_string += stuff` to combine strings since `join()` takes linear time rather than O(_n_^2) time.



## From Google's Python Style Guide

See https://google.github.io/styleguide/pyguide.html

* Use mypy, pylint, and/or PyCharm inspections.
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
* A function must have a docstring, unless it's short, obvious, and not externally visible. Explain what you need to know to call it.
* Classes should have docstrings.
* Don't put code at the top level that you don't want to run when the file is imported, unit tested, or pydoc'd.



## Other Recommendations

* Prefer `format_string.format(...)` over `format_string % ...` because [printf-style % formatting has](https://docs.python.org/3/library/stdtypes.html?highlight=sprintf#printf-style-string-formatting) "a variety of quirks that lead to a number of common errors (such as failing to display tuples and dictionaries correctly)". `.format()` is a bit more readable and it offers some nice additional features. (See: More info on [the 4 ways to format strings in Python](https://dbader.org/blog/python-string-formatting).)

# Performance Tips

[Summarized from sources like [PythonSpeed](https://wiki.python.org/moin/PythonSpeed).]

* Testing membership in a set or a dict is very fast, `O(1)`, unlike a list, tuple, or array.
* Sorting a list using a _sort key_ is faster than using a _comparison function._
* Mapping a function over a list, or using a list comprehension or generator comprehension, should be faster than a `for` loop since it pushes the loop work into compiled C code.
* Local variables are faster to access than global variables, builtins, and attribute lookups.
* Iterators are generally more memory-friendly and scalable than list operations, e.g. `xrange()`, `itertools.imap()`, `dict.iteritems` vs. `range()`, `map()`, `dict.items()`.
* Core building blocks are coded in optimized C, including the builtin datatypes (lists, tuples, sets, and dictionaries) and extension modules like `array`, `itertools`, and `collections.deque`.
* Builtin functions run faster than hand-built equivalents, e.g. `map(operator.add, v1, v2)` is faster than `map(lambda x, y: x+y, v1, v2)`.
* For queue applications using `pop(0)` or `insert(0, v)`, `collections.deque()` offers superior `O(1)` performance over a list because it avoids the `O(n)` step of rebuilding a list for each insertion or deletion.
* Chained comparisons like `x < y < z` are faster and hopefully more readable than `x < y and y < z`.
* Threading can improve the response time in applications that would otherwise waste time waiting for I/O.
