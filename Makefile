clean:
	rm -rf build/
	rm -rf dist/
	rm -rf .eggs/
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +
test: 
	pytest

publish: package
	twine upload dist/*

package: clean 
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist
