#!/bin/bash

rsync -avP -e 'ssh -i /home/dstrauss/.ssh/ec2key.rsa' sgeadmin@ec2-50-16-29-100.compute-1.amazonaws.com:subsurface/src/contrastX .