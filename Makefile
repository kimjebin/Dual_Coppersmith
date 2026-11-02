.PHONY: check tables smoke
check:
	python scripts/verify_submission.py
	python -m unittest discover -s tests -v
tables:
	python scripts/make_tables.py
smoke:
	python scripts/smoke.py --work /tmp/dual-coppersmith-smoke
