export LD_LIBRARY_PATH=/stor9000/apps/users/NWSUAF/2018055070/sf/glibc_2.14/lib
## NCBI v5 nr database addition: taxonomy DB
~/sf/diamond makedb --threads 20 --in nr_191101.fa --db nr.dmnd --masking 0 --taxonmap prot.accession2taxid --taxonnodes ./taxdmp/nodes.dmp --taxonnames ./taxdmp/names.dmp
