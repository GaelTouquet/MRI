export GEOMETRY=$PWD
export PYTHONPATH=$GEOMETRY/..:$PYTHONPATH
export PATH=$GEOMETRY:$PATH

# Install all dependencies listed in the requirements.txt file
pip install -r requirements.txt
