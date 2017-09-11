#This file contains the code behind Example 4.3
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#The files "SphericalProjections.gap" and "ProjectiveReduction.gap"
#should be available through the same source as this file.


#Create the Mathieu group M_11 acting on 12 points
G:=TransitiveGroup(12,272);;


#Find a stabilizer for the action on pairs of poinst
H:=Stabilizer(G,[1,2],OnPairs);


#Examine the constituents of the permutation character to find the
#position of the one with degree 10
chi:=PermutationCharacter(G,H);;
const:=ConstituentsOfCharacter(chi);;
List(const,psi->DegreeOfCharacter(psi));

#It's the second one


#Compute the spherical projections
Read("SphericalProjections.gap");
e:=SphericalProjections(G,H);;


#Combine the projections for the trivial character and the one of degree 10
p:=e[1]+e[2];;


#Projectively reduce the result
Read("ProjectiveReduction.gap");
p:=ProjectiveReduction(p);;


#Rescale to make it a projection again
p:=2*p;