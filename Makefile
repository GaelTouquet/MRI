install:
	pip3 install -r requirements.txt
format:
	autopep8 --in-place *.py
example:
	python3 ./temporal_resolution.py
