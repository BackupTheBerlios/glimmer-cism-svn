#!/bin/bash
#
# Magnus Hagdorn
# April 2010
#
# script for backing up glimmer-cism subversion repo


# check if a keychain is available for this computer
if [ -f ${HOME}/.keychain/$(hostname)-sh ]; then
   source ${HOME}/.keychain/$(hostname)-sh
fi

mkdir -p ${HOME}/gc-backup

/usr/bin/ssh svn.berlios.de "svnadmin -q dump /svnroot/repos/glimmer-cism | bzip2" > ${HOME}/gc-backup/glimmer-cism.dump.bz2
