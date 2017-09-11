#This file contains the code behind Example 4.4
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#The files "SphericalProjections.gap" and "ProjectiveReduction.gap"
#should be available through the same source as this file.


#Construct the groups involved
sl25:=SL(2,5);;
L:=SylowSubgroup(sl25,2);;
gl25:=GL(2,5);;
G0:=Normalizer(gl25,L);;
G:=SemidirectProduct(G0,GF(5)^2);;


#Note: GAP stores G as a matrix group, where
#a pair (M,x), with M in G0 and x in GF(5)^2,
#is represented by a matrix like [[M,0],[x,1]]
#
#To get the affine action of G on a point
#y in GF(5)^2, we treat the latter as [y[1],y[2],1]
#instead.


#Get the stabilizer of a pair of points in GF(5)^2
H:=Stabilizer(G,[ Z(5)^0*[1,0,1], Z(5)^0*[0,1,1] ],OnPairs);;


#To verify that G acts doubly transitive, it suffices to
#check that Size(G)/Size(H) = 25*24, the left hand side
#being the size of the orbit of the pair stabilized by H,
#and the right hand side being the total number of pairs
#of distinct points in GF(5)^2.
Size(G)/Size(H);


#Create the groups K0 and K, and verify that the latter is
#regular for the action of G on X.

K0:=Normalizer(sl25,L);;
StructureDescription(K0);

K:=SemidirectProduct(K0,GF(5)^2);;

#It has trivial intersection with H:
Size(Intersection(Elements(H),Elements(K)));

#And G=HK:
Size(G) = Size(H)*Size(K);


#Find the constituents of the permutation character that have
#degree 2 and are complex conjugates of each other
chi:=PermutationCharacter(G,H);;
const:=ConstituentsOfCharacter(chi);

#We want the second and the third.

#Verify that they occur with multiplicity one
ScalarProduct(chi,const[2]);
ScalarProduct(chi,const[3]);


#Compute the spherical projections
Read("SphericalProjections.gap");
e:=SphericalProjections(G,H);;


#Make the 2x6
Read("ProjectiveReduction.gap");
p26:=ProjectiveReduction(e[2]);;

#Make the 4x12
p412:=ProjectiveReduction(e[2]+e[3]);;

#Rescale to make them projections again
p26:=300*2/6*p26;

p412:=150*4/12*p412;