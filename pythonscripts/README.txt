---------------------------------------------------------------
Usage help:: 
	python viewRun.py --help 
	python viewAssembledSum.py --help

--------------------------------------------------------------
Getting PYTHON to work on SLAC's Unix computers:

To obtain the correct PYTHON libraries, you could run the following shell scripts provided by Andy Salnikov:

test -f /reg/g/psdm/etc/ana_env.sh && . /reg/g/psdm/etc/ana_env.sh

You could append this to your $HOME/.bashrc file to always run these shell scripts when you login your shell. When you do this for the first time, you will need either run "source $HOME/.bashrc" (or its equivalent) or restart a new 
shell to enforce these changes (logout, then login again). 

---------------------------------------------------------------

