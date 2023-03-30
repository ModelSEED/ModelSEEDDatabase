#!/bin/bash
/home/seaver/Software/Java_11/bin/javac -cp "/home/seaver/Software/MarvinBeans/lib/*:." -Xlint:deprecation $1
cp $1 ../../../Structures/Mol/
