#!/bin/bash

VERSION=0.1
mkdir -p ${HOME}/.vmdplugins/dbtools${VERSION}
rsync -av --exclude="*.git" . ${HOME}/.vmdplugins/dbtools${VERSION}
