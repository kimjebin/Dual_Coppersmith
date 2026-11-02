"""Use the original successive-minimum proof with resumable sphere partitions."""
import certify
import partition_search
certify.exact_engine._search=partition_search.search
certify.EXACT_SEARCH_LABEL='partitioned fplll 5.5.0 error-bounded complete enumeration and exact integer lifts'
if __name__=='__main__': certify.main()
