export MRI=$PWD
export PYTHONPATH=$MRI/..:$PYTHONPATH
export PATH=$MRI:$PATH

# Install all dependencies listed in the requirements.txt file
pip3 install -r requirements.txt
