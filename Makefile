.PHONY: test

test:
	nosetests -v -s deicode --with-coverage --cover-package=deicode
