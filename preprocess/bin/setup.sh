echo 'Setting up...'

read -p "If you have not created and activated a python3 virtual environment, \
now is the time to do so. Continue? [y/n]" yn
case $yn in
    [Yy] ) continue;;
    [Nn] ) exit;;
    * ) echo "Please answer y or n.";;
esac
echo
echo 'Installing Python3 modules...'
pip install -r requirements.txt --quiet
echo 'Successfully installed Python3 modules.'
echo
echo 'Creating data directory and subdirectories...'
mkdir data
mkdir data/raw
mkdir data/intermediate
mkdir data/preprocessed
mkdir data/postprocessed
echo 'Created data directory and subdirectories.'
echo 'Done.'

