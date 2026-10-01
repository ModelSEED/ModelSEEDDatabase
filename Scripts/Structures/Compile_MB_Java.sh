#!/bin/bash
/home/seaver/Software/Java_11/bin/javac -cp "/home/seaver/Software/MarvinBeans/lib/*:." -Xlint:deprecation $1
file=`basename $1 .java`
cp ${file}.class ../../../Structures/Mol/
