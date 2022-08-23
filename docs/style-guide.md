# Style Guide

_These guidelines are condensed from the PEP8 and Google style guides here, along with our decisions. -- Jerry_

* See [PEP8](https://www.python.org/dev/peps/pep-0008/)
* See [Google's Python style guide](https://google.github.io/styleguide/pyguide.html)
* See our [March 12, 2018 meeting slides](https://docs.google.com/presentation/d/1pf6GQmMwbUoeASmNk1bmjYZ-PteJ0tSFH0v6P6vd5eE/edit#slide=id.g313d94100c_0_55) with notes on the tradeoffs between alternatives.

## Note: Python 2 and 3 compatibility
The project has transitioned to Python 3.

We have not yet removed Python 2 compatibility elements such as
`from __future__ import absolute_import, division, print_function`
and uses of the `six` compatibility library.


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

* Use pylint and/or PyCharm inspections.
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


# Type hints

Python type hints can help catch bugs, improve code documentation (what to pass in this argument?),
and aid development tools such as PyCharm's code completion and refactoring. They can also provide input to
foreign language bridges like Jython to Java, DB query mapping, and RPC parameter marshalling. Reports
are that developers love it when their coworkers add type hints.

Python will remain dynamically typed. Type hints will not become mandatory even by convention.

The key advice is to use them to tease apart bytes vs. unicode since that's the core gotcha in the
transition to Python 3.

This is a good place to apply the 80/20 rule. Roughly speaking, do the easy 20% of the work to find 80% of
the type bugs. Complicated typing cases can get very complicated in programming languages. So mainly put
type hints on functions and occasionally on variables when the type checker can't figure it out, e.g. mypy
can't infer the element type when you create an empty collection, but pytype can figure it out by looking
ahead to the code that adds elements. In a difficult case, you can punt by using type `Any` which can
assign to or from anything.

See [these Google Slides](https://docs.google.com/presentation/d/1xwVHjpQsRTGLdaIpPDNjKPQFtp7mBNfn7oTWK4Kct30/edit?usp=drivesdk)
for our plan overview for migrating to Python 3 and adding type hints.

Until we finish the transition to Python 3, use Python 2+3 compatible type hint syntax like the
following examples.


## Examples

```python
def emphasize(message, suffix=''):
  # type: (str, str) -> str
  """Construct an emphatic message. Docstring comes after the function type hint."""
  return message + '!' + suffix
```

In the example above, `suffix` is an optional argument to the caller but it's always a `str` to the
type checker, so we declare it as an ordinary `str` here. `Optional[str]` means "str or None" as in
the following example:

```python
from typing import Optional

def title(name, suffix=None):
  # type: (str, Optional[str]) -> str
  """Construct an emphatic message."""
  return name + ('' if suffix is None else ' ' + suffix)
```

```python
from typing import Dict

def headline(
        text,           # type: str
        width=80,       # type: int
        fill_char="-",  # type: str
    ):                  # type: (...) -> str
    for x, y in points:  # type: float, float
        something = g(x, y)  # type: ignore
    return "This is a way to handle lots of arguments."

def f(x, *args, **kwargs):
    # type: (str, *float, **str) -> Dict[str]
    """For *args, give the type of an arg value, not the tuple.
    For **kwargs, give the type of an arg value, not the dict.
    """
    return kwargs

class C(object):
    def __init__(self, p):
        # type: (float) -> None
        """Don't give a type for self or cls. __init__() return type is None."""
        self.p = p

    def circumference(self, radius):
        # type: (float) -> float  # account for all args _except_ `self` or `cls`
        version = (2, 7, 14)  # type: tuple
        return 2 * math.pi * radius
```


## The Typing module

To retrofit static typing tools to existing versions of Python 2, the designers added
* the type hint comment syntax,
* [the "typing" module](https://docs.python.org/3/library/typing.html) on [PyPI](https://pypi.org)
  containing types like `List` and `Dict` with uppercase names so they don't have to patch the
  existing classes `list` and `dict`,
* `.pyi` "stub files" to add type definitions onto existing Python libraries and native libraries
  (see below for numpy type stubs),
* the [Typeshed](https://github.com/python/typeshed) repository for "stub" files
  (it's bundled with PyCharm, mypy, and pytype)

You can use class names like `dict`, `list`, and `tuple`, but that won't specify the content types.

* `None`: the value None (this is easier than writing `NoneType`)
* `Any`: any type; assigns to or from anything
* `Dict[str, float]`: a dict of str to float, that is floats keyed with strings
* `Mapping[str, float]`: an abstract mapping; accepts a `dict`, `OrderedDict`, or other
* `List[int]`: a list of int values
* `Tuple[int, int, int]`: a tuple containing 3 integers
* `Optional[str]`: a str or None
* `Callable[[str, List[int]], int]`: a function or anything else you can call given a str and a
list of ints, returning an int
* `Union[str, int]`: accepts either a str or an int
* `Sequence[float]`: any sequence of floats that supports `.__len__()` and `.__getitem__()`
* `Iterable[float]`: any iterable collection of floats; accepts a `list`, `tuple`, 
* `Text`: a unicode string in both Python 2 and 3
* `AnyStr`: any kind of text or byte string but doesn't allow different kinds of strings to mix


## Tools

[PyCharm](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html) checks types interactively
while you edit. You can also run a particular kind of inspection (such as "Type checker")
or all inspection types on a single file or any source directory. See
[Python Type Checking (Guide)](https://realpython.com/python-type-checking/).

[mypy](https://github.com/python/mypy/) is a batch program to check types. We run it in
Continuous Integration. Code changes must run cleanly to pass the `ecoli-pull-request`
and the `ecoli-small` CI tests. To run it:

    mypy

[PyAnnotate](https://github.com/dropbox/pyannotate) (Py2.7) and
[MonkeyType](https://github.com/Instagram/MonkeyType) (Py3) will observe types at runtime and
write them out as stub files or proposed type hints in the source code.


## Types for Numpy

There are experimental type stubs in the numpy repo [numpy-stubs](https://github.com/numpy/numpy-stubs)
that define types for `dtype` and `ndarray`. It's not fancy but it does catch some mistakes and it
improves PyCharm autocompletion. Hopefully the numpy team will improve these stubs, but numpy is more
flexible with types than the type system is unlikely to handle.

With this stub file, you can write type hints like `np.ndarray`, `np.ndarray[int]`, and `np.ndarray[Any]`.
It doesn't have a way to express array shape so the shape still goes into a docstring.

```python
import numpy as np

def f(a):
    # type: (np.ndarray[float]) -> np.ndarray[int]
    return np.asarray(a, dtype=int)
```


The wcEcoli project includes numpy-stubs.

To install more stub files:
1. Copy them into a `stubs/` directory in the project.
2. Mark the `stubs/` directory as a source root in PyCharm by choosing **Mark Directory as | Sources Root**
from the directory's context menu.


## Tips

* When a method overrides a superclass method, it inherits the type hints. You needn't
repeat them but if you do, they must match.
* Call `reveal_type(xxx)` to ask the type checker to print its inferences.
* Escape hatches: `Any`, `cast()`, `# type: ignore`, and `.pyi` stub files (esp. for C extensions).
* Gradually add types to code, one file at a time, starting with the most heavily used code.
* Run the type checker in development and in Continuous Integration. Fix its warnings to defend progress.
We can tell mypy to disallow untyped functions in particular modules.
* The config setting `check_untyped_defs = True` will also check the contents of functions that don't
have types hints. It might find actual bugs.
* Lambdas don't support type hints, so `def`ine a function when you want typing.
* To punt on types for a difficult functions, just leave out the type annotations or add a
  `@no_type_check` decorator. That treats it as having the most general type possible, i.e. each arg
  type is `Any` (except a method's `self` or `cls` arg) and the result is `Any`. `@no_type_check`
  might also disable type inferences and checking inside the function.


## Terminology

* A _type_ is for type checking, variable annotations, and function annotations.
  A _class_ is a runtime thing. Every _class_ acts like a _type_, and there are additional types like
  `Union[str, int]` which you can't use as classes.
* Type A is a _subtype_ of B if it has a subset of the values (classification) **and** a superset of the
methods. A value of a subtype can act like a value from its supertypes.
* The type `From` _is consistent with_ (assignable to a variable of) type `To` if
  * `From` is a subtype of `To`, _or_
  * `From` or `To` is `Any`.
* _Type hints_ just enable tools: docs, type checkers, IDE code completion, etc.
* The type checker only complains about _inconsistent_ types. It does nothing at runtime.
* Type _annotations_ do the same job as type hint comments, but only in Python 3.
  These annotations are available at runtime via the `__annotations__` attribute.
* _Gradual typing_ means adding type hints to existing code, freely mixing code
  with and without have type hints.
* The type checker can _infer the types_ of local variables, `@property` methods, and more.
* _Nominal typing_ is based on class **names**. Python types are mostly nominal.
* _Structural typing_ is based on structure such as the duck-type protocol `Typing.Sized`
  which means "it has a `.__len__()` method".
* Type information on values is _erased_ at runtime.


## Covariant and Contravariant types

This goes beyond the 80/20 rule but you might need to know about it. See
[The Ultimate Guide to Python Type Checking](https://realpython.com/python-type-checking/) for a
good tutorial that includes this stuff.

Some generic types (_sources_ like Tuple[t1, t2] and FrozenSet[t]) are _covariant_, some are _contravariant_
(_sinks_ like Callable[[t1], …]), and some are _invariant_ (like List[t]).

   * `bool` is a subtype (and a subclass) of `int`.
   * `Tuple[bool]` is a subtype of `Tuple[int]`. It's _covariant_ with the element type.
   * `Callable[[bool], None]` is a supertype of `Callable[[int], None]`. It's _contravariant_.
   * `List[bool]` is _invariant_ with `List[int]`. It might look like a subtype but `mylist.append(10)`
   would make it not a `List[bool]`.
   * `T = TypeVar('T', covariant=True)` declares covariance for a generic via its type variables.
   * `T = TypeVar('T', contravariant=True)` types are invariant by default.


## References

[PEP 484: Type Hints](https://www.python.org/dev/peps/pep-0484/). Parts of this PEP are worth reading.

[The "typing" pip](https://docs.python.org/3/library/typing.html).

[Python type checking tools](https://python-type-checking.readthedocs.io/en/latest/tools.html)
   * mypy [overview](http://mypy-lang.org/), [docs](http://mypy.readthedocs.org),
   [blog](http://mypy-lang.blogspot.com), [github](github.com/python/mypy).
   * PyCharm
     * It can help [add type hints for you](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html):
     Press Opt-Return (or click the "intention" lightbulb icon) then "Add type hint for ..."
     * [Type Hinting in PyCharm](https://www.jetbrains.com/help/pycharm/type-hinting-in-product.html)
     shows how to use PyCharm features to quickly add type hints, validate them, and eventually convert
     comment-based type hints to Py3.6 variable annotations once we leave Py2 behind.
   * [pytype](https://github.com/google/pytype) Google's static analyzer checks types, attribute names, and more.
   * [Pyre](github.com/facebook/pyre-check) Facebook's static checker is faster on large code bases than mypy
   * [pyannotate](github.com/dropbox/pyannotate) runtime type info collector
   * [MonkeyType](github.com/Instagram/MonkeyType) runtime type info collector


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
