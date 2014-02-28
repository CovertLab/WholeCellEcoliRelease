**Framework**
 * determine whether to split unique instances away from MoleculeCounts
  * agreed to separate the classes
 * remove old unique instances code from class
  * Nick is working on this
 * rename MoleculeCounts class and other references
  * Nick is working on this
 * refactor class to be more in line with current needs (after all the other work is done)

**Knowledge base**
 * types of attributes (simple numbers, object ID references, states)
 * adding unique instance attrs to KB (can be fake for now)
 * list-of-arrays instantiation from KB
  * currently using a placeholder

**Accessing**
 * attribute type accessors
  * accesors for non-basic types are deferred until we have a need
 * indexing rules (adding, removing, extending allocated space)
  * complete for everything other than merging partitions

**Querying (partition setup and requests)**
 * query format
  * finished
 * passing queries
  * finished
 * evaluating queries
  * finished

**Paritioning**
 * trivial partitioning (only one process requests)
  * Nick is working on a process
 * deferred partitioning (handling by other states)
  * needs a new state
 * request collision handling
  * needs a good algorithm

**Disk**
 * table creation
  * finished
 * table appending
  * finished
 * table loading
  * finished

**Misc**
 * build actual Chromosome State
 * Replication Process
 * update Transcription/Translation (clean up logic, add comments, change to TUs instead of gene/protein pairs)
 * Transcripts State (better name?)
 * Write summary of state model behavior
