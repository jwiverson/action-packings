#This file contains the code behind Example 3.6
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#It requires the package "FinInG", which at the time of writing was
#available at http://cage.ugent.be/fining/
#
#The file "SphericalProjections.gap" should be available through the
#the same source as this file.


LoadPackage("fining");

#Create the affine space AG(3,2) and make its affine linear group
as:=AffineSpace(3,2);;
G:=AffineGroup(as);;


#We want the action of G on lines in AG(3,2), so we compute
#the stabilizer of a random line
H:=Stabilizer(G,Random(Lines(as)));;


#Compute the spherical projections
Read("SphericalProjections.gap");
e:=SphericalProjections(G,H);;


#Use trace to check the ranks of the projections,
#and find the one with rank 7
List(e,p->Trace(p));

#It's the third one

p:=e[3];