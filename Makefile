clean:
	find . -name "*.cPickle" -exec rm -fr {} \;
	find . -name "*.pyc" -exec rm -rf {} \;
