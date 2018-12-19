## Comparing ways to assemble byte strings into a 2D NumPy array

There are non-obvious performance impacts of various ways for `TableReader.readColumn2D()` to decompress and unpack a sequence of data blocks into a 2D NumPy array.

Code outline:
* Read the data blocks.
* Decompress them into byte strings.
* Unpack each byte string into a 2-D NumPy array with one or many rows (depending on how many rows were packed into the block).
* Use the given `indices` (if any) to select subcolumns of the arrays.
* Stack the arrays into a taller 2-D array.

Strategies and lessons learned:
1. `join()` the byte strings together, then convert them to a NumPy array via `np.frombuffer(str).reshape(...)`. That lets native code loops do most of the work. But some callers require the resulting array to be writable, and `np.frombuffer(str)` returns a read-only array onto that read-only byte `str`. Selecting subcolumns of the output array will get a writable copy, otherwise the code had better copy it.
2. Unpack each block into an array, then call `np.vstack()` to combine them into one (writable) array. However it makes another copy of the data in RAM alongside the compressed blocks and the unpacked arrays. For large arrays this turns out to be slow.
3. Allocate the output array via `np.zeros()`, then copy each unpacked array into the correct rows. By keeping each decompressed byte string in memory only while unpacking and copying that data to the result array, this can save lots of RAM and some time. However the code can't allocate the result array until it knows the total number of rows. Alternatives:
   1. Decompress and unpack the blocks, then calculate the total number of rows.
   2. The header says how many rows are in each block, although the last block might not be full. So collect the compressed blocks, then decompress the last block to count its rows, then allocate the result array, decompressing only one block at a time to save most of the working RAM and ≈9% of the read time.
   3. Have the writer backpatch the total number of rows into the file header. This saves ≈14% of the read time at the cost of writer complexity and time, and if the writer process exits before this completion step, the reader can't read the files (unless it implements the above strategy as a fallback).
   4. Have the writer write the total number of rows into a table attribute. This makes column files depend on the table's attributes (unless the reader implements the fallback strategy).
   5. Have the writer write a chunk containing the total row count after the last data block. The reader would collect up compressed blocks, then read the count chunk and allocate the result array. This doesn't really save anything except it replaces the special case for the last block with a special/error case for a missing count chunk.
4. Use `np.empty()`, but measurements show it's no faster than `np.zeros()`.

See:
* tablereader.py, tablewriter.py
* measure_bulk_reader.py
* measure_zlib.py


### Measuring ways to combine the data

The following code times the ways of combining the decompressed data blocks. `measure()` is called with a list of bytestrings from a simulation output, along with the target array dtype. Two of these ways turned out to return read-only arrays.

These measurements are on a 2018 MacBook Pro with Core i9 CPU.

```
from time import clock
def timeit(f): c = clock; start = c(); f(); end = c(); return end - start

def f1(bs, dt):  # read-only result
  return np.frombuffer(b''.join(bs), dt).reshape(len(bs), -1)

def f2(bs, dt):
  return np.vstack([np.frombuffer(b, dt) for b in bs])

def f3(bs, dt):
  result = np.empty((len(bs), len(bs[0]) // dt.itemsize), dt)
  for i, b in enumerate(bs):
    result[i, :] = np.frombuffer(b, dt)
  return result

def f4(bs, dt):  ###
  result = np.zeros((len(bs), len(bs[0]) // dt.itemsize), dt)
  for i, b in enumerate(bs):
    result[i, :] = np.frombuffer(b, dt)
  return result

def f5(bs, dt):  ###
  return np.frombuffer(b''.join(bs), dt).reshape(len(bs), -1).copy()

def f6(bs, dt):  # read-only result
  return np.frombuffer(b''.join(bs), dt).reshape(len(bs), -1)[:, :]

def measure(bytelist, dtype_):
  m = [
    timeit(lambda: f1(bytelist, dtype_)),
    timeit(lambda: f2(bytelist, dtype_)),
    timeit(lambda: f3(bytelist, dtype_)),
    timeit(lambda: f4(bytelist, dtype_)),
    timeit(lambda: f5(bytelist, dtype_)),
    timeit(lambda: f6(bytelist, dtype_))
  ]
  print('{:9.6f}, {:9.6f}, {:9.6f}, {:9.6f}, {:9.6f}, {:9.6f}'.format(*m))

BulkMolecules/atpRequested:
shape = (2901, 12), nbytes = 2901 * 96 = 278496
   f1 --      f2 #4      f3 #2      f4 #3     f5 #1!     f6 --
 0.000086,  0.005782,  0.002312,  0.002343,  0.000058,  0.000039,
 0.000092,  0.006948,  0.002388,  0.002383,  0.000059,  0.000039,
 0.000068,  0.005936,  0.002313,  0.002306,  0.000059,  0.000042,
 0.000075,  0.006022,  0.002374,  0.002385,  0.000058,  0.000040,
 0.000071,  0.005912,  0.002438,  0.002433,  0.000061,  0.000044
 --------   --------   --------   --------   --------   --------
 0.000078   0.006120   0.002365   0.002370   0.000059   0.000041

BulkMolecules/counts:   measure(cb, counts.dtype)
shape = (2901, 38575), nbytes = 2901 * 308600 = 895248600
   f1 --      f2 #2      f3 #3      f4 #1      f5 #4      f6 --
 0.575963,  0.579574,  0.591050,  0.565404,  1.203193,  0.564041,
 0.580189,  0.581060,  0.585343,  0.581092,  1.173112,  0.553914,
 0.578382,  0.598036,  0.595418,  0.586048,  1.175091,  0.582533,
 0.575067,  0.575303,  0.592868,  0.576914,  1.169940,  0.562360,
 0.569895,  0.580466,  0.582237,  0.594195,  1.191714,  0.569332,
 --------   --------   --------   --------   --------   --------
 0.5758992  0.5828878  0.5893832  0.5807306  1.18261    0.566436

Summary:
* f2, f3, f4 are about the same time for large 'counts'.
* f3, f4 are about the same, and faster than f2 for 'atpRequested'.
```

**Note:** The run times are very different for the modest size `BulkMolecules/atpRequested` column vs. the large `BulkMolecules/counts` column.


### Q. Why is the BulkMolecules/counts so much slower? Swapping?

```
# bc = byte strings from BulkMolecules/counts
# dt = its dtype

>>> timeit(lambda: f4(bc, dt))
0.6148289999999994
>>> timeit(lambda: f4(bc, dt))
0.6127919999999998
>>> timeit(lambda: f4(bc[:len(bc)/4], dt))*4
0.6046960000000006
>>> timeit(lambda: f4(bc[:len(bc)//5], dt))*5
0.6249150000000014
>>> timeit(lambda: f4(bc[:len(bc)//10], dt))*10
0.6294100000000036
>>> timeit(lambda: f4(bc[:len(bc)//100], dt))*100
0.635700000000039

>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.6024800000000141
>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.33692000000002054
>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.332080000000019
>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.3365200000000357
>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.3446399999999983
>>> timeit(lambda: f4(bc[:len(bc)//80], dt))*80
0.6611200000000395
>>> timeit(lambda: f4(bc[:len(bc)//80], dt))*80
0.17424000000005435
>>> timeit(lambda: f4(bc[:len(bc)//80], dt))*80
0.16760000000004993
>>> timeit(lambda: f4(bc[:len(bc)//80], dt))*80
0.1695200000000341
>>> timeit(lambda: f4(bc[:len(bc)//80], dt))*80
0.1512800000000425
>>> timeit(lambda: f4(bc[:len(bc)//40], dt))*40
0.34019999999998163

>>> ll = len(bc)//40*40
>>> timeit(lambda: f4(bc[:ll], dt))
0.6114519999999999

>>> ll = ll/40; ll
73
>>> timeit(lambda: f4(bc[:ll], dt)) * 40
0.33472000000003277
>>> timeit(lambda: f4(bc[:ll], dt)) * 40
0.37916000000002725
>>> timeit(lambda: f4(bc[:ll], dt)) * 40
0.3117600000000209
```

A. It doesn't speed up much more than proportionally for small memory cases, so apparently not.


### Timing readColumn2D variations

Several variations of `readColumn2D()` were line-profiled side-by-side on BulkMolecules/counts, BulkMolecules/atpRequested, UniqueMoleculeCounts/time, and some other table columns from an E. coli simulation.

* Collecting up the compressed bytes rather than the decompressed and unpacked arrays saved `readColumn2D()` ≈9% of the overall time relative to the `vstack()` approach on `counts`, perhaps due to caching.
* Simulating the reader knowing the total number of rows up front -- as if the writer back-patched it into the column header -- the reader could allocate the result array then unpack and copy in each block as soon as it was read. That saved maybe ≈15% of the time relative to the `vstack()` approach, but there is a lot of measurement noise.
* A slight variation on the combining loop releases each compressed block ASAP. This profiled about the same, maybe a hair slower, but definitely within measurement noise:

```
while entry_blocks:  # instead of `for raw in entry_blocks:`
  entries = decomp(entry_blocks.pop(0))
  ...
```
Another idea is to use `np.memmap` to lazily load uncompressed data.  Early work (see Issue #221) suggests that this may lead to cryptic performance issues.
