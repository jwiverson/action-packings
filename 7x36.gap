#This file contains the Gram matrix of the 7x36 frame constructed
#in Example 4.2 of the paper "Optimal Line Packings From Finite
#Group Actions", by Joseph W. Iverson, John Jasper, and Dustin G.
#Mixon.


g:=[ [ 1, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -3/7, -3/7, 3/7, -3/7, 3/7, 3/7, -1/7, -3/7, -3/7, -3/7, 3/7, 3/7, 1/7, 1/7, 1/7, 1/7, 1/7, 
      -3/7, -3/7, -3/7, 3/7, -3/7, -3/7, 3/7, -3/7, -3/7, -3/7 ], [ -1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, -3/7, 3/7, 
      3/7, -3/7, -3/7, 3/7, 3/7, 3/7, -3/7, -3/7, -3/7, -3/7, 3/7, 3/7, 3/7, -3/7, 3/7, -3/7, 3/7, -3/7, 3/7 ], 
  [ -1/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, 3/7, -3/7, 3/7, 3/7, -1/7, -3/7, -3/7, -3/7, 3/7, -3/7, 3/7, -1/7, -3/7, -3/7, -1/7, -3/7, 
      3/7, 3/7, 3/7, -1/7, -3/7, 3/7, -1/7, -3/7, -3/7, 1/7, 1/7, 3/7 ], [ -1/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, 3/7, -1/7, -3/7, -3/7, 3/7, -3/7, 3/7, 
      1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 3/7, -3/7, 3/7, 3/7, -3/7, 3/7, 3/7, -3/7, -3/7, -3/7, 3/7, 3/7, -3/7, 3/7, 3/7 ], 
  [ -1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, -3/7, -3/7, 3/7, -3/7, 3/7, 3/7, -1/7, 3/7, -3/7, 3/7, -1/7, -3/7, -3/7, 3/7, -3/7, -3/7, 
      -1/7, 3/7, 3/7, -3/7, 3/7, -1/7, 3/7, -1/7, 3/7, -1/7, -3/7, 1/7 ], [ -1/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, 3/7, 3/7, -1/7, -3/7, -3/7, 3/7, -3/7, 
      3/7, -1/7, -3/7, 3/7, -3/7, 3/7, -1/7, 3/7, 3/7, -3/7, 3/7, 1/7, 1/7, 1/7, 1/7, 3/7, 3/7, -3/7, -3/7, 3/7, -3/7 ], 
  [ -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, -3/7, 3/7, 3/7, -1/7, -3/7, -3/7, 3/7, -3/7, -3/7, 3/7, 3/7, 3/7, -1/7, -3/7, -3/7, 3/7, 
      -3/7, -1/7, -3/7, 3/7, -1/7, -3/7, 3/7, -3/7, -1/7, 3/7, -1/7, -1/7 ], [ -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, -3/7, 3/7, -3/7, 3/7, 3/7, -1/7, 
      -3/7, 3/7, 3/7, -1/7, -3/7, 3/7, -3/7, -3/7, 3/7, -1/7, 3/7, -3/7, -1/7, -3/7, 3/7, -3/7, 1/7, 1/7, 1/7, 3/7, 3/7, -3/7 ], 
  [ -1/7, 1/7, 3/7, 3/7, -3/7, 3/7, -3/7, -3/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 3/7, -3/7, 3/7, -3/7, 3/7, 1/7, 1/7, 1/7, 
      1/7, 1/7, 3/7, 3/7, -3/7, 3/7, -3/7, 3/7, -3/7, -3/7, 3/7, 3/7 ], [ -3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 3/7, 3/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, 1/7, 
      1/7, 1/7, 1/7, 1/7, 1/7, -3/7, 3/7, 3/7, -3/7, -3/7, -3/7, 3/7, 3/7, -3/7, 3/7, 3/7, -3/7, 3/7, 3/7, -3/7 ], 
  [ -3/7, 1/7, 3/7, -3/7, 3/7, -1/7, 3/7, -3/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, -3/7, -1/7, 3/7, 3/7, -3/7, -3/7, -1/7, -3/7, -3/7, 
      -3/7, 3/7, 1/7, 1/7, 1/7, 1/7, 3/7, -3/7, -3/7, 3/7, -3/7, 3/7 ], [ 3/7, 1/7, 3/7, -3/7, -3/7, -3/7, -1/7, 3/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, 
      -3/7, 3/7, -3/7, -3/7, 3/7, -1/7, -3/7, 3/7, -3/7, 3/7, -1/7, -3/7, -3/7, -1/7, 3/7, -3/7, -3/7, -1/7, 3/7, -1/7, -1/7 ], 
  [ -3/7, 1/7, -1/7, 3/7, 3/7, -3/7, -3/7, 3/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, 3/7, 3/7, 3/7, -3/7, -1/7, -3/7, 3/7, -1/7, -3/7, 
      3/7, -3/7, 3/7, -1/7, 3/7, -3/7, -1/7, 3/7, 3/7, 1/7, 1/7, 3/7 ], [ 3/7, 1/7, -3/7, -3/7, 3/7, 3/7, -3/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, 3/7, 
      -3/7, -1/7, -3/7, -3/7, 3/7, 3/7, 3/7, -1/7, -3/7, 3/7, -1/7, -3/7, 3/7, 3/7, 1/7, 1/7, 1/7, -3/7, -3/7, -3/7 ], 
  [ 3/7, 1/7, -3/7, 3/7, -1/7, -3/7, 3/7, -3/7, -1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, -3/7, -3/7, 3/7, -1/7, 3/7, 3/7, 3/7, -3/7, 3/7, 
      -1/7, -3/7, -3/7, 3/7, -3/7, -1/7, -3/7, -1/7, 3/7, -1/7, -3/7, 1/7 ], [ -1/7, -3/7, -3/7, 1/7, 3/7, 3/7, -3/7, 3/7, -1/7, 1/7, -3/7, -3/7, 3/7, 3/7, -3/7, 
      1, -1/7, -1/7, -1/7, -1/7, -1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 3/7, -3/7, 3/7, -3/7, 3/7, 3/7, 3/7, -3/7, 3/7, -3/7 ], 
  [ -3/7, 3/7, 3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 3/7, 1/7, -1/7, 3/7, 3/7, -3/7, -3/7, -1/7, 1, -1/7, -1/7, -1/7, -1/7, -1/7, 3/7, -3/7, 
      3/7, -3/7, 1/7, 1/7, 1/7, 1/7, -3/7, 3/7, -3/7, 3/7, 3/7, 3/7 ], [ -3/7, 3/7, -3/7, 1/7, 3/7, -3/7, 3/7, -1/7, -3/7, 1/7, 3/7, -3/7, 3/7, -1/7, 3/7, -1/7, 
      -1/7, 1, -1/7, -1/7, -1/7, 3/7, -3/7, -1/7, -3/7, -3/7, -1/7, 3/7, 3/7, -3/7, 1/7, 1/7, 1/7, 3/7, -3/7, 3/7 ], 
  [ -3/7, -3/7, 3/7, 1/7, -1/7, 3/7, 3/7, -3/7, 3/7, 1/7, 3/7, -3/7, -3/7, -3/7, -1/7, -1/7, -1/7, -1/7, 1, -1/7, -1/7, -3/7, -3/7, 3/7, 
      -1/7, 3/7, 3/7, 3/7, -3/7, -1/7, 3/7, -1/7, -3/7, -1/7, 3/7, 1/7 ], [ 3/7, -3/7, -1/7, 1/7, -3/7, -3/7, 3/7, 3/7, -3/7, 1/7, -3/7, 3/7, -1/7, -3/7, 3/7, 
      -1/7, -1/7, -1/7, -1/7, 1, -1/7, -3/7, -1/7, 3/7, 3/7, -3/7, -3/7, -1/7, -3/7, -3/7, -1/7, -3/7, 3/7, 1/7, 1/7, -3/7 ], 
  [ 3/7, 3/7, -3/7, 1/7, -3/7, 3/7, -1/7, -3/7, 3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 3/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, 3/7, 3/7, 3/7, 
      -3/7, -1/7, -3/7, 3/7, -1/7, 3/7, -3/7, 3/7, -1/7, -3/7, -1/7, -1/7 ], [ 1/7, 3/7, -3/7, 3/7, 3/7, -1/7, -3/7, -3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 3/7, 3/7, 
      1/7, -1/7, 3/7, -3/7, -3/7, 3/7, 1, -1/7, -1/7, -1/7, -1/7, 1/7, 1/7, 1/7, 1/7, -3/7, 3/7, 3/7, -3/7, -3/7, 3/7 ], 
  [ 1/7, 3/7, -1/7, -3/7, -3/7, 3/7, -3/7, 3/7, 1/7, 3/7, -3/7, 3/7, -1/7, 3/7, -3/7, 1/7, 3/7, -3/7, -3/7, -1/7, 3/7, -1/7, 1, -1/7, 
      -1/7, -1/7, -3/7, -1/7, 3/7, 3/7, -1/7, 3/7, -3/7, 1/7, 1/7, -3/7 ], [ 1/7, -3/7, -3/7, 3/7, -3/7, 3/7, 3/7, -1/7, 1/7, 3/7, -3/7, -3/7, -3/7, -1/7, 3/7, 
      1/7, -3/7, -1/7, 3/7, 3/7, 3/7, -1/7, -1/7, 1, -1/7, -1/7, -1/7, 3/7, -3/7, -3/7, 1/7, 1/7, 1/7, -3/7, 3/7, -3/7 ], 
  [ 1/7, -3/7, 3/7, 3/7, -1/7, -3/7, -3/7, 3/7, 1/7, -3/7, -3/7, 3/7, 3/7, -3/7, -1/7, 1/7, 3/7, -3/7, -1/7, 3/7, -3/7, -1/7, -1/7, -1/7, 
      1, -1/7, 3/7, -3/7, -3/7, -1/7, -3/7, -1/7, 3/7, -1/7, 3/7, 1/7 ], [ 1/7, -3/7, 3/7, -3/7, 3/7, 3/7, -1/7, -3/7, 1/7, -3/7, 3/7, -1/7, -3/7, 3/7, -3/7, 
      1/7, -3/7, -3/7, 3/7, -3/7, -1/7, -1/7, -1/7, -1/7, -1/7, 1, 3/7, -3/7, -1/7, 3/7, 3/7, -3/7, -1/7, -3/7, -1/7, -1/7 ], 
  [ -3/7, -3/7, 3/7, 3/7, 3/7, 1/7, -3/7, -1/7, 3/7, -3/7, 1/7, -3/7, 3/7, -1/7, -3/7, 3/7, 1/7, -1/7, 3/7, -3/7, -3/7, 1/7, -3/7, -1/7, 
      3/7, 3/7, 1, -1/7, -1/7, -1/7, 1/7, 1/7, 1/7, -3/7, 3/7, 3/7 ], [ -3/7, 3/7, -1/7, 3/7, -3/7, 1/7, 3/7, -3/7, 3/7, 3/7, 1/7, -3/7, -1/7, -3/7, 3/7, -3/7, 
      1/7, 3/7, 3/7, -1/7, 3/7, 1/7, -1/7, 3/7, -3/7, -3/7, -1/7, 1, -1/7, -1/7, -1/7, 3/7, -3/7, 1/7, 1/7, 3/7 ], 
  [ -3/7, 3/7, -3/7, -3/7, 3/7, 1/7, -1/7, 3/7, -3/7, 3/7, 1/7, -1/7, 3/7, 3/7, -3/7, 3/7, 1/7, 3/7, -3/7, -3/7, -1/7, 1/7, 3/7, 
      -3/7, -3/7, -1/7, -1/7, -1/7, 1, -1/7, 3/7, 3/7, -1/7, 3/7, -1/7, -1/7 ], [ 3/7, 3/7, 3/7, -3/7, -1/7, 1/7, -3/7, -3/7, 3/7, -3/7, 1/7, 3/7, -3/7, 3/7, 
      -1/7, -3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 1/7, 3/7, -3/7, -1/7, 3/7, -1/7, -1/7, -1/7, 1, -3/7, -1/7, -3/7, -1/7, -3/7, 1/7 ], 
  [ -3/7, -3/7, -1/7, -3/7, 3/7, 3/7, 3/7, 1/7, -3/7, 3/7, 3/7, -3/7, -1/7, 1/7, -3/7, 3/7, -3/7, 1/7, 3/7, -1/7, -3/7, -3/7, -1/7, 1/7, 
      -3/7, 3/7, 1/7, -1/7, 3/7, -3/7, 1, -1/7, -1/7, 1/7, 1/7, -3/7 ], [ -3/7, 3/7, -3/7, 3/7, -1/7, 3/7, -3/7, 1/7, 3/7, 3/7, -3/7, -3/7, 3/7, 1/7, -1/7, 3/7, 
      3/7, 1/7, -1/7, -3/7, 3/7, 3/7, 3/7, 1/7, -1/7, -3/7, 1/7, 3/7, 3/7, -1/7, -1/7, 1, -1/7, -1/7, 3/7, 1/7 ], 
  [ 3/7, -3/7, -3/7, 3/7, 3/7, -3/7, -1/7, 1/7, -3/7, -3/7, -3/7, -1/7, 3/7, 1/7, 3/7, 3/7, -3/7, 1/7, -3/7, 3/7, -1/7, 3/7, -3/7, 
      1/7, 3/7, -1/7, 1/7, -3/7, -1/7, -3/7, -1/7, -1/7, 1, -3/7, -1/7, -1/7 ], [ -3/7, 3/7, 1/7, -3/7, -1/7, -3/7, 3/7, 3/7, -3/7, 3/7, 3/7, 3/7, 1/7, -3/7, 
      -1/7, -3/7, 3/7, 3/7, -1/7, 1/7, -3/7, -3/7, 1/7, -3/7, -1/7, -3/7, -3/7, 1/7, 3/7, -1/7, 1/7, -1/7, -3/7, 1, -1/7, 1/7 ], 
  [ -3/7, -3/7, 1/7, 3/7, -3/7, 3/7, -1/7, 3/7, 3/7, 3/7, -3/7, -1/7, 1/7, -3/7, -3/7, 3/7, 3/7, -3/7, 3/7, 1/7, -1/7, -3/7, 1/7, 3/7, 
      3/7, -1/7, 3/7, 1/7, -1/7, -3/7, 1/7, 3/7, -1/7, -1/7, 1, -1/7 ], [ -3/7, 3/7, 3/7, 3/7, 1/7, -3/7, -1/7, -3/7, 3/7, -3/7, 3/7, -1/7, 3/7, -3/7, 1/7, -3/7, 
      3/7, 3/7, 1/7, -3/7, -1/7, 3/7, -3/7, -3/7, 1/7, -1/7, 3/7, 3/7, -1/7, 1/7, -3/7, 1/7, -1/7, 1/7, -1/7, 1 ] ];;