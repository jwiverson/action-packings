#This file contains the code behind Example 4.2
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#The files "SphericalProjections.gap" and "ProjectiveReduction.gap"
#should be available through the the same source as this file.


G:=SL(2,8);;

#Create a new method, OnPairsOfLines, that allows G to act accordingly
OnPairsOfLines:=function(pair,g)
	return [OnLines(pair[1],g), OnLines(pair[2],g)];
end;;


#Create a base point for this action
x0:=[ Z(8)^0*[1,0], Z(8)^0*[0,1] ];;


#Find the stabilizer of x0
H:=Stabilizer(G,x0,OnPairsOfLines);;


#Make the permutation character and analyze its constituents
chi:=PermutationCharacter(G,H);;
const:=ConstituentsOfCharacter(chi);;

#To see that exactly one of the constituents is real of degree 7, type
#
#Browse(const);
#
#It's const[2]

#Verify that it occurs with multiplicity one
ScalarProduct(chi,const[2]);


#Compute the spherical projections
Read("SphericalProjections.gap");
e:=SphericalProjections(G,H);;


#Projectively reduce e[2]
Read("ProjectiveReduction.gap");
p:=ProjectiveReduction(e[2]);;


#Rescale to make it a projection again
p:=2*p;


