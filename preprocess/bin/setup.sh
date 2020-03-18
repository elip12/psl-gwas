echo 'Setting up...'

echo "If you have not created and activated a python3 virtual environment,"
echo "now is the time to do so."
read -p "Continue? [y/n] " yn
case $yn in
    [Yy] ) ;;
    [Nn] ) exit;;
    * ) echo "Please answer y or n.";;
esac
echo 'Installing Python3 modules...'
pip install -r requirements.txt --quiet
echo 'Successfully installed Python3 modules.'
echo 'Creating data directory and subdirectories...'
mkdir data
mkdir data/raw
mkdir data/intermediate
mkdir data/preprocessed
mkdir data/postprocessed
echo 'Successfully created data directory and subdirectories.'
echo 'Done.'

