"""Use proved error-bounded enumeration; retain the original two-pass logic."""
from shared_bootstrap import *
from shared import exact_minima
import certify
from fast_search import search

exact_minima._search = search
certify.EXACT_SEARCH_LABEL = 'fplll error-bounded exhaustive enumeration with exact integer lifts'
certify.main()
