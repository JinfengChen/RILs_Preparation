virtualenv env
source env/bin/activate
export PYTHONPATH=`pwd`/env/lib/python2.7/site-packages
pip install --upgrade pip
pip install "numpy"
pip install "scipy"
pip install "matplotlib"
pip install "pandas"
pip install "Bio"
pip install "Biopython"
#run script here
deactivate env
