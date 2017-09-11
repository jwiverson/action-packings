#This file contains the code behind Example 4.1
#of the paper "Optimal Line Packings From Finite Group Actions",
#by Joseph W. Iverson, John Jasper, and Dustin G. Mixon.
#
#The file "SphericalFunctions.gap" should be available through
#the same source as this file.

#Create the generators of the group K
t:=[[0,1],[1,0]];;
m:=[[1,0],[0,-1]];;
i:=[[1,0],[0,1]];;

T1:=KroneckerProduct(KroneckerProduct(t,i),i);;
T2:=KroneckerProduct(KroneckerProduct(i,t),i);;
T3:=KroneckerProduct(KroneckerProduct(i,i),t);;
M1:=KroneckerProduct(KroneckerProduct(m,i),i);;
M2:=KroneckerProduct(KroneckerProduct(i,m),i);;
M3:=KroneckerProduct(KroneckerProduct(i,i),m);;
iI:=E(4)*IdentityMat(8);;

K:=Group(T1,T2,T3,M1,M2,M3,iI);;


#Create the normalizing matrices U and V, and
#the group H that they generate

U:=E(8)/Sqrt(2)*[[0,  0,  -1,  0,  E(4),  0,  0,  0],
[0,  0,  -E(4),  0,  1,  0,  0,  0],
[0,  0,  0,  E(4),  0,  1,  0,  0],
[0,  0,  0,  1,  0,  E(4),  0,  0],
[-1,  0,  0,  0,  0,  0,  E(4),  0],
[E(4),  0,  0,  0,  0,  0,  -1,  0],
[0,  E(4),  0,  0,  0,  0,  0,  1],
[0,  -1,  0,  0,  0,  0,  0,  -E(4)]];;

V:=E(8)/Sqrt(2)*[[0,  0,  0,  0,  E(4),  -1,  0,  0],
[0,  0,  0,  0,  -E(4),  -1,  0,  0],
[E(4),  1,  0,  0,  0,  0,  0,  0],
[-E(4),  1,  0,  0,  0,  0,  0,  0],
[0,  0,  E(4),  -1,  0,  0,  0,  0],
[0,  0,  -E(4),  -1,  0,  0,  0,  0],
[0,  0,  0,  0,  0,  0,  E(4),  1],
[0,  0,  0,  0,  0,  0,  -E(4),  1]];;

H:=Group(U,V);;


#Verify that H is isomorphic to PSU(3,3)
StructureDescription(H);


#Verify that H normalizes K
IsNormal(H,K);


#Verify that H intersects K trivially
Size(Intersection(Elements(H),Elements(K)));


#Make the group G=SemidirectProduct(H,K)
G:=Group(T1,T2,T3,M1,M2,M3,iI,U,V);;


#To verify that (G,H) is a Gelfand pair, read 
#Prof. Akihide Hanaki's file "association_scheme.gap", 
#which at time of writing was available at 
# http://math.shinshu-u.ac.jp/~hanaki/as/
#
#Read("association_scheme.gap");
#
#Then type
#
#IsCommutativeScheme(SchurianScheme(G, H));


#Make the table of spherical functions
Read("SphericalFunctions.gap");

tbl:=SphericalFunctions(G,H);
